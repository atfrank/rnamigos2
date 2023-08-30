import pickle
import sys
import os
from collections import Counter
import itertools

import networkx as nx
from tqdm import tqdm
import dgl
import numpy as np
import torch
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader, Subset
from dgl.dataloading import GraphDataLoader
import pandas as pd

from rnaglib.utils import NODE_FEATURE_MAP 


class DockingDataset(Dataset):
    def __init__(self,
                 annotated_path,
                 edge_types=None,
                 nuc_types=None,
                 shuffle=False,
                 target='dock',
                 seed=0,
                 debug=False
                 ):
        """
            Setup for data loader.

            Arguments:
                annotated_path (str): path to annotated graphs (see `annotator.py`).
                get_sim_mat (bool): whether to compute a node similarity matrix (deault=True).
                nucs (bool): whether to include nucleotide ID in node (default=False).
        """
        print(f">>> fetching data from {annotated_path}")
        self.path = annotated_path
        self.all_graphs = sorted(os.listdir(annotated_path))
        if debug:
            self.all_graphs = self.all_graphs[:100]
        if seed:
            print(f">>> shuffling with random seed {seed}")
            np.random.seed(seed)
        if shuffle:
            np.random.shuffle(self.all_graphs)
        #build edge map
        self.edge_map = {e: i for i, e in enumerate(sorted(edge_types))}
        self.num_edge_types = len(self.edge_map)
        self.target = target

        self.n = len(self.all_graphs)


    def __len__(self):
        return self.n


    def __getitem__(self, idx):
        """
            Returns one training item at index `idx`.
        """
        _, graph, _, ring, fp_nat, fp, inter_score, inter_score_trans, score_native_ligand, label_native_lig, label_1std, label_2std, label_thr_min30, label_thr_min17, label_thr_min12, label_thr_min8, label_thr_0, sample_type, is_native  = pickle.load(open(os.path.join(self.path, self.all_graphs[idx]), 'rb'))

        #adding the self edges
        # graph.add_edges_from([(n, n, {'label': 'X'}) for n in graph.nodes()])
        graph = nx.to_directed(graph)
        one_hot = {edge: torch.tensor(self.edge_map[label.upper()]) for edge, label in
                   (nx.get_edge_attributes(graph, 'label')).items()}
        nx.set_edge_attributes(graph, name='edge_type', values=one_hot)

        node_attrs = None
        one_hot_nucs  = {node: NODE_FEATURE_MAP['nt_code'].encode(label) for node, label in
                (nx.get_node_attributes(graph, 'nt')).items()}

        nx.set_node_attributes(graph, name='nt_features', values=one_hot_nucs)

        #g_dgl = dgl.DGLGraph()
        g_dgl = dgl.from_networkx(nx_graph=graph, edge_attrs=['edge_type'], node_attrs=['nt_features'])
        g_dgl.title = self.all_graphs[idx]

        if self.target == 'fp':
            target = fp
        if self.target == 'dock':
            target = inter_score_trans

        return g_dgl, fp_nat, torch.tensor(target, dtype=torch.float), [idx]

class Loader():
    def __init__(self,
                 dataset,
                 batch_size=128,
                 num_workers=20,
                 shuffle=False,
                 seed=0,
                 debug=False
                 ):
        """
        Wrapper class to call with all arguments and that returns appropriate data_loaders
        :param pocket_path:
        :param ligand_path:
        :param batch_size:
        :param num_workers:
        :param augment_flips: perform numpy flips
        :param ram: store whole thing in RAM
        :param siamese: for the batch siamese technique
        :param full_siamese for the true siamese one
        """
        self.dataset = dataset
        self.batch_size = batch_size
        self.num_workers = num_workers
        self.seed = seed
        self.shuffle = shuffle


    def get_data_strat(self, k_fold=0):
        n = len(self.dataset)
        indices = list(range(n))
        labels = []
        for idx, item in enumerate(self.dataset):
            labels.append(item[2])

        # collate_block = collate_wrapper()

        if k_fold > 1:
            from sklearn.model_selection import StratifiedKFold
            kf = StratifiedKFold(n_splits=k_fold)
            #from sklearn.model_selection import KFold
            #kf = KFold(n_splits=k_fold)
            for train_indices, test_indices in kf.split(np.array(indices), np.array(labels)):
                train_set = Subset(self.dataset, train_indices)
                test_set = Subset(self.dataset, test_indices)

                train_loader = GraphDataLoader(dataset=train_set, shuffle=True, batch_size=self.batch_size,
                                          num_workers=self.num_workers, collate_fn=None)
                test_loader = GraphDataLoader(dataset=test_set, shuffle=True, batch_size=self.batch_size,
                                         num_workers=self.num_workers, collate_fn=None)

                yield train_loader, test_loader

        else:
            from sklearn.model_selection import train_test_split

            train_indices, test_indices, train_indices_labels, test_indices_labels = train_test_split(np.array(indices), np.array(labels), test_size=0.2,
                                                           train_size=0.8, random_state=None,
                                                           shuffle=True, stratify=np.array(labels))

            train_set = Subset(self.dataset, train_indices)
            test_set = Subset(self.dataset, test_indices)

            train_loader = GraphDataLoader(dataset=train_set, shuffle=True, batch_size=self.batch_size,
                                      num_workers=self.num_workers, collate_fn=None)
            test_loader = GraphDataLoader(dataset=test_set, shuffle=True, batch_size=self.batch_size,
                                 num_workers=self.num_workers, collate_fn=None)

            # return train_loader, valid_loader, test_loader
            yield train_loader, test_loader

    def get_data(self):
        train_size = int(0.8 * len(self.dataset))
        test_size = len(self.dataset) - train_size
        train_set, test_set = torch.utils.data.random_split(self.dataset, [train_size, test_size])

        train_loader = GraphDataLoader(dataset=train_set, shuffle=self.shuffle, batch_size=self.batch_size,
                                  num_workers=self.num_workers, collate_fn=None)
        test_loader = GraphDataLoader(dataset=test_set, shuffle=self.shuffle, batch_size=self.batch_size,
                             num_workers=self.num_workers, collate_fn=None)

        return train_loader, test_loader


def describe_dataset(annotated_path='../data/annotated/pockets_docking_annotated'):
    rdock_scores = pd.DataFrame(columns=['POCKET_ID', 'LABEL_NAT', 'LABEL_1STD', 'LABEL_2STD', 'TOTAL'])
    path = annotated_path
    all_graphs = sorted(os.listdir(annotated_path))
    for g in tqdm(all_graphs):
        #p = pickle.load(open(os.path.join(path, g), 'rb'))
        _, graph, _, ring, fp_nat, _, fp, total_score, label_nat, label_1std, label_2std = pickle.load(open(os.path.join(path, g), 'rb'))
        rdock_scores = rdock_scores.append({'POCKET_ID': g, 
            'LABEL_NAT': str(label_nat), 
            'LABEL_1STD': str(label_1std), 
            'LABEL_2STD': str(label_2std), 
            'TOTAL': str(total_score)}, ignore_index=True)

    pickle.dump(rdock_scores, open('dataset_labels_and_score.p', 'wb'))
    rdock_scores.to_csv('dataset_labels_and_score.csv')



class InferenceLoader(Loader):
    def __init__(self,
                 annotated_path,
                 batch_size=5,
                 num_workers=20):
        super().__init__(
            annotated_path=annotated_path,
            batch_size=batch_size,
            num_workers=num_workers)
        self.dataset.all_graphs = sorted(os.listdir(annotated_path))

    def get_data(self):
        collate_block = collate_wrapper()
        train_loader = GraphDataLoader(dataset=self.dataset,
                                  shuffle=False,
                                  batch_size=self.batch_size,
                                  num_workers=self.num_workers,
                                  collate_fn=None)
        return train_loader


if __name__ == '__main__':
    loader = Loader(shuffle=False,seed=99, batch_size=1, num_workers=1)
    data = loader.get_data(k_fold=1)
    pass
