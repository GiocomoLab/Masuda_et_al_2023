Run the jupyter notebook in an environment which contains all necessary packages. Eg in Conda:
````
$ conda create -n ketamine python=3.9
$ conda activate ketamine
(ketamine) conda install -c anaconda jupyter
(ketamine) pip install spectral_connectivity spikeinterface h5py phylib
(ketamine) jupyter notebook
````
