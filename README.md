# MiR_Combo_Targeting
Companion code for the paper: "Network potential identifies therapeutic miRNA cocktails in Ewing sarcoma"
You can find the pre-print for this analysis here.https://www.biorxiv.org/content/10.1101/854695v1
If you have any questions regarding how to run this pipeline on your device/ use it for your research, 
feel free to shoot me an email at davis.weaver@case.edu. 

## Here is the general flow: 

Most of the heavy lifting is done by an R package that I wrote entitled "disruptr". 
You can install disruptr using the command: 
```
devtools::install_github("https://github.com/DavisWeaver/disruptr")
```

### You can skip this step if you like! 

The main network potential pipeline can then be run by sourcing the file `Main_EWS_MiR.R`. 
This will likely require a computing cluster (if you comment out the `compute_null` command it can run on a laptop in a few hours). 
Sample batch files are provided for submitting to a linux cluster.

### Start here if you are short on time

I have included the output of our main computational pipeline as a system file in the disruptr package. 
To run the second part of our analysis (the miRNA optimization bit), simply open `miR_Analysis_Ewing.R` and hit run. 
After that, the `miR_Figs_Ewing.R` script that generates all of our figures and tables should run mostly without a hitch. 
If you would like to reproduce all of our figures - you will need to get access to some data that we are not authorized to distribute (namely the transcriptomic data from St. Jude Children's Research Hospital)


## Dependencies:

Package dependencies are at the top of each script. 
`miR_Figs_Ewing.R` will break into a million pieces if you try and run it before `miR_Analysis_Ewing.R`. 
One of our figures makes use of the data provided by this paper: https://pubmed.ncbi.nlm.nih.gov/31978347/. We did not include it in this repository because it is too large and github didn't like it but it is freely available for download from the original authors.

# Questions

Feel free to reach out to me at davis.weaver@case.edu if you have any questions! If you would like to contribute / learn more about _in-silico_ repression techniques, head on over to https://github.com/DavisWeaver/disruptr 

