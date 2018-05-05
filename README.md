# A Shiny server for Bioinformatics apps
A web service accessible through a web browser (from anywhere), that holds a collection of **shiny** apps to perform routine bioinformatics analyses.

## shinotate
A shiny server app to deliver processed and annotated transcriptome data stored in [Trinotate](https://trinotate.github.io/) db

### Main features
* Provide basic summarized information (sequences, annotation, expression)
* Provide basic plotting (expression) 

### Implementation
#### Hosting
* Shiny server hosted on a cloud-service (AWS, NECTAR, DigitalOcean, etc.) - see [details](http://deanattali.com/2015/05/09/setup-rstudio-shiny-server-digital-ocean/)
* Mount USC's main servers using `sshfs` - see [details](https://www.digitalocean.com/community/tutorials/how-to-use-sshfs-to-mount-remote-file-systems-over-ssh) 

<<<<<<< HEAD
### Homepage
* Select from available transcriptomes 
* Provide basic information

### Data presentation
#### Metadata
* Experimental design
* Analysis methods (Rmarkdown)
#### Transcript information
=======
#### Transcript information
* Select from available transcriptomes  
>>>>>>> 25bceb4dbc99dc179b97039e50395ee1e2363f9d
* Allow filtering based on keyword/transcript/gene  
* Download selected sequences or annotation table

#### Expression
* Show selected transcript expression (in a table and a graph)
* Allow consolidating multiple transcripts under one gene name/symbol 

## CFX-qpcR
A shiny server to analyse qPCR results from Biorad CFX platform.

### Main features
* Select genes/targets
* Warn about outliers/extreme variability among replicates
* Display calibration curves plots
* Save linear regression, R<sup>2</sup> and calculated *Eff* in a table available for download
