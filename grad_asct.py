import torch
import torch.nn as nn
import numpy as np
import argparse
import pandas as pd
from argparse import ArgumentParser
import importlib.util

## NEED TO CHANGE THESE FILE PATHS
spec = importlib.util.spec_from_file_location("happynet.models", "/home/jar268/happynet_grph/happynet_graph/happynet/models.py")
f = importlib.util.module_from_spec(spec)
spec.loader.exec_module(f)

spec2 = importlib.util.spec_from_file_location("utils.data_utils", "/home/jar268/happynet_grph/happynet_graph/happynet/utils/data_utils.py")
f2 = importlib.util.module_from_spec(spec2)
spec2.loader.exec_module(f2)



#TO DO: GET THE REST OF THE MODEL TYPES RUNNING
IMPLEMENTED_MODELS = [#'lstm',
                      #'lstm_aux',
                      'seq2seq_aux',
                      'cnn',
                      #'cnn_ae',
                      #'cnn_ext'
                      ]

####
def s2s_enc(seq_csv,loaded_model):
	test_seqs = pd.read_csv(seq_csv)
	# tbatch = test_seqs['CDR3']
	tbatch = test_seqs.iloc[:,0]
	cnt = 0
	for seq in tbatch:
		tbatch[cnt] = [f2.seq_to_inds(seq)]
		cnt +=1

	emb = torch.cat(tuple(torch.tensor(i) for i in tbatch),0)
	embedded_batch = loaded_model.forward(emb)[1]
	print(embedded_batch.shape)
	return(embedded_batch)

def cnn_enc(seq_csv,loaded_model):
	test_seqs = pd.read_csv(seq_csv)
	tbatch = test_seqs.iloc[:,0]
	cnt = 0
	for seq in tbatch:
		tbatch[cnt] = [f2.seq_to_inds(seq)]
		cnt +=1
	emb = torch.cat(tuple(torch.tensor(i) for i in tbatch),0)
	return(emb)

def s2s_dec(embedding,num_layers=2,num_directions=2,hidden_size=50,seq_len=20):
	hidden_state = embedding.reshape(-1, num_layers * num_directions, hidden_size)

	# decoder expecting shape (num_layers * num_directions, batch, hidden_size
	hidden_state = hidden_state.transpose(0,1)

	#option 2
	cell_state = torch.zeros(hidden_state.shape)


	enc_out = torch.randn(embedding.shape[0], seq_len, hidden_size)
	enc_out = enc_out.type_as(hidden_state)

	dec_out = loaded_model.decode(enc_out, (hidden_state, cell_state))

	#translate from one-hot to a.a. seq
	m = nn.Softmax(2) #compute logsoftmax along seq dimension 
	soft = m(dec_out)
	maxx = torch.argmax(soft,dim=2)

	out = [f2.inds_to_seq(i.tolist()) for i in maxx] #convert indices to amino acids
	cnt = 0
	for i in out:
		out[cnt] = [''.join(i)]
		cnt+=1
	return(out)



## this currently performs an optimization loop (k=100) for a single sequence/latent space coordinate.
def optimize(embeddings,k=100,lr=1e-6):

	#embeddings should be shape (batch size,200)

	for emb in range(embeddings.shape[0]): #optimization loop for each sequence embedding
		with torch.enable_grad():
			embeddings = embeddings.type(torch.float)
			print(embeddings[emb,])
			input_data = embeddings[emb,].requires_grad_(True)
			for step in range(k):

				if cl_args.model== 'cnn':
					out = model(input_data)[0][1]
				if cl_args.model == 'seq2seq_aux':
					out = model(input_data)


				if not input_data.requires_grad:
					raise ValueError("input not gradient-enabled")
				if not out.requires_grad:
					raise ValueError("output not gradient-enabled")


				criterion = nn.MSELoss(reduction='sum')
				loss = criterion(out,input_data)
				grad = torch.autograd.grad(loss,input_data)

				input_data2 = input_data.clone()
				input_data2 += grad[0]*lr

				input_data = input_data2 #this is because pytorch doesn't allow in-place operations on leaf variables w/ grad=T

				if step == 0 or step == k-1:
					print('loss, rd{}: '.format(step) + str(loss.item()))
			# print(input_data)

		embeddings[emb] = input_data

	return(embeddings)




if __name__ == "__main__":

	parser = ArgumentParser(description="pass the type of model, the run location (where embeddings are located, if needed), the statedict location, and the seed sequences")
	parser.add_argument('--model', default='seq2seq_aux', type=str)
	parser.add_argument('--run_location', default='/home/jar268/happynet_grph/happynet_graph/seq2seq_aux_giffordseq2seq_aux/gifford/2020-10-18-14-39-48/wandb/run-20201018_183948-hdd7a746/happynet_project/hdd7a746/checkpoints/epoch=76.ckpt', type=str)
	parser.add_argument('--statedict', default='/home/jar268/happynet_grph/happynet_graph/optim/enrichment_model', type=str)
	parser.add_argument('--seed', default='seed_seqs.csv', type=str)

	cl_args = parser.parse_args()
	assert cl_args.model in IMPLEMENTED_MODELS, 'model not implemented'


	## PRE-TRAINED MODEL WEIGHTS LOADING

	m = f.str2model(cl_args.model)

	if cl_args.model == 'seq2seq_aux':
		##TO DO: PULL HPARAMS FROM SAVED MODEL
		hparamss = argparse.Namespace(lr=0.0001,embedding_dim=50,hidden_dim=50,layers=2,probs=0.2,bidirectional=True,reg_ramp=True,alpha_val=0.5)
		seq2seq = m(hparamss)
		loaded_model = seq2seq.load_from_checkpoint(cl_args.run_location)
		#####
		print('pre-trained model loaded')
		## enrichment_model is state_dict created from MLP pre-trained on latent space
		model = nn.Sequential(nn.Linear(200,32),nn.ReLU(),nn.Linear(32,32),nn.ReLU(),nn.Linear(32,1))
		model.load_state_dict(torch.load(cl_args.statedict))
		print('regressor loaded')
		embedded_batch = s2s_enc(cl_args.seed,loaded_model)


	if cl_args.model == 'cnn':
		hparamss = argparse.Namespace(lr=0.0001,embedding_dim=25,hidden_dim=50,layers=2,probs=0.2,bidirectional=True,reg_ramp=True,alpha_val=0.5,input_dim=22,kernel_size=4,seq_len=237)
		cnn = m(hparamss)
		loaded_model = cnn.load_from_checkpoint(cl_args.run_location)
		#####
		print('pre-trained model loaded')
		model = loaded_model
		print('no regressor')
		embedded_batch = cnn_enc(cl_args.seed,loaded_model)


	
	opt = optimize(embedded_batch)
	out = s2s_dec(opt)

	pd.DataFrame(out).to_csv('optimized_sequences.csv')







