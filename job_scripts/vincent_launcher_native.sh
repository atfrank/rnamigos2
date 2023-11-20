#!/bin/bash

CMDARRAY=()

for split in {0..9};
do
	for dim in 16 64;
    do
    python_cmd="python experiments/train.py
                      model.encoder.hidden_dim=16
                      model.encoder.hidden_dim=${dim}
                      model.decoder.in_dim=$((dim+32))
                      model.decoder.out_dim=1
                      model.decoder.activation=sigmoid
					            model.dropout=0.3
                      train.target=is_native
                      train.loss=bce
                      train.num_epochs=1000
                      train.early_stop=100
                      train.learning_rate=1e-3
                      train.rnamigos1_split=${split}
                      name=rnamigos2_native_dim${dim}_${split}
                      device=cpu
                      train.num_workers=0"
    python_cmd=$(echo $python_cmd) # to replace newlines
    CMDARRAY+=("$python_cmd")
    done
done

for dim in 16 64;
  do
  python_cmd="python experiments/train.py
                  model.encoder.hidden_dim=16
                  model.encoder.hidden_dim=${dim}
                  model.decoder.in_dim=$((dim+32))
                  model.decoder.out_dim=1
                  model.decoder.activation=sigmoid
                  model.dropout=0.3
                  train.target=is_native
                  train.loss=bce
                  train.num_epochs=1000
                  train.early_stop=100
                  train.learning_rate=1e-3
                  name=native_chembl_dim${dim}
                  device=cpu
                  train.num_workers=0"
done
python_cmd=$(echo $python_cmd) # to replace newlines
CMDARRAY+=("$python_cmd")

echo "${CMDARRAY[13]}"

#N_JOBS=50
#for ((i = 0; i < ${#CMDARRAY[@]}; i++))
#do
#    eval "${CMDARRAY[$i]}" &
#    [ $(( i+1 % N_JOBS )) -eq 0 ]  && wait
#done


