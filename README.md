# FTM
## Files
- Codes: FTM.zip
- Examples: examples.zip
- Datasets: BJDATA.zip EURDATA.zip EURDATA.z01 SYNDATA.zip

## Environment
- IDE: VS2017
- Operation System: win10 x64
***

## Parameters of Execution Program
- i : File name of input temporal graph (default: example-temporal-graph.txt)
- f : File name of output (default: example-output.txt)
- k : File name of frequency threshold (default:example-k.txt, if not exists k=10), the program can test different frequency thresholds for one time (untested)

```
file content format:
     10
     20
     30
```

- r : Algorithm id (default:1)
	
	- 4:static algorithm FTM
	- 9:incremental algorithm DFTM

- o : Output level (default:0)
	
&ensp;&ensp; three levels:
		
	- 0:only output motif number, the running time and memory use
		
	- 1:except for those outputs mentioned above, output the edges number of every motif
		
	- 2:except for those outputs mentioned above, output the detailed information of motif edges(nodes and labels)

- g : Index id (default:1)

	- 1: eL table
	- 5: IC tree

- l : Limit the start time and end time of input data (default:doesn't limit) format:-l:0,500

e.g.: if the interval of input temporal graph is [1,100], -l:1,50 means you only use the interval [1,50] of temporal graph when testing the algorithm  (used in incremental algorithm)

- n : Limit the end time of input data when snapshot increasing (use - l at the same time) (default:doesn't limit) 

e.g.: -l:0,500 -n:1600 means that the interval of temporal graph is [0,500] before snapshots increase, and that the interval of temporal graph is [0,1600] after snapshots increase (used in incremental algorithm)

- e : File name of fixed labels of motif edges (default: do not fix labels)

```
 file content format (means motifs with label types set = {1,2}, all edges have label type in {1,2}, one line for one label type):
     1
     2				
```
- m : Method to save the set TF (default:2)
	
	- 1 : maintain the temporal motifs in intervals (the output are in the time order in one row of TI-Table; need O(T^2) space; we use this method for sorted outputs)
	- 2 : maintain the temporal motifs in rows (the output are not in the time order in one row of TI-Table; does not need O(T^2) space) 

**Examples for static algorithm FTM**

1. FTNMotifs.exe -i:temporalgraph.txt -k:k10.txt -r:4 -g:1 -f:output.txt -o:1  

 (use EL table, only output the number of motifs and edges, time, memory)

2. FTNMotifs.exe -i:temporalgraph.txt -k:k10.txt -r:4 -g:5 -f:output.txt -o:2 

  (use IC trees, output nodes of edges, labels, the number of motifs and edges, time, memory) 

**Example for incremental algorithm DFTM**

3. FTNMotifs.exe -i:temporalgraph.txt -k:k10.txt -r:9 -g:1 -f:output.txt -o:0 -l:0,500 -n:1000    

(use EL table, only output the number of motifs, time, memory, interval of temporal graph is from [0,500] to [0,1000]) 
***

## Datasets
- BJD<font size = 8>ATA</font>: A real-life dataset records road traffic conditions in Beijing. There are three traffic conditions (i.e. 2: congested, 1: slow, -1: fast), and are updated every 5 minutes. 
- EURD<font size = 8>ATA</font>: A real-life dataset for a renewable European electric power system. Each edge represents a transmission line with only static properties, and each node represents a merge point of transmission lines with a dynamic property of the hourly energy demand. 
- LAND<font size = 8>ATA</font>: A real-life dataset for real-life computer event dataset from the [Los Alamos National Laboratory enterprise network](https://csr.lanl.gov/data/2017/#citing). Each node represents a computer, and each edge represents the communication of two computers with a dynamic property of the average  number of bytes sent.
- SYND<font size = 8>ATA</font>: Synthetic datasets produced by the synthetic data generator.