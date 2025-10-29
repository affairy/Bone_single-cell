#软骨
import os,sys
os.getcwd()
os.listdir(os.getcwd())

import loompy as lp;
import numpy as np;
import scanpy as sc;
x=sc.read_csv("scenic.data.csv")
row_attrs = {"Gene": np.array(x.var_names),}
col_attrs = {"CellID": np.array(x.obs_names)}
lp.create("scenic.loom",x.X.transpose(),row_attrs,col_attrs)


###滑膜成纤维
import os,sys
os.getcwd()
os.listdir(os.getcwd())

import loompy as lp;
import numpy as np;
import scanpy as sc;
x=sc.read_csv("scenic.data.csv")
row_attrs = {"Gene": np.array(x.var_names),}
col_attrs = {"CellID": np.array(x.obs_names)}
lp.create("scenic.loom",x.X.transpose(),row_attrs,col_attrs)