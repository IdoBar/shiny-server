# A Shiny server for Bioinformatics apps
A web service accessible through a web browser (from anywhere), that holds a collection of **shiny** apps to perform routine bioinformatics analyses.

## shinotate
A shiny server app to deliver processed and annotated transcriptome data stored in [Trinotate](https://trinotate.github.io/) db

### Main features
* Provide basic summarized information (sequences, annotation, expression)  
* Provide basic plotting (expression)  

### Implementation
#### Hosting
* Shiny server hosted on a cloud-service (AWS, **NECTAR**, DigitalOcean, etc.) - see [details](http://deanattali.com/2015/05/09/setup-rstudio-shiny-server-digital-ocean/)  
* Enable user authentication and secure connection (TLS/SSL) - either implemented through `nginx` (see this [guide](https://www.r-bloggers.com/add-authentication-to-shiny-server-with-nginx/)), or can be hosted by a third-party service, such as [Auth0](https://auth0.com/blog/adding-authentication-to-shiny-server/)  
* Mount analysis servers to the cloud-server using `sshfs` - see [details](https://www.digitalocean.com/community/tutorials/how-to-use-sshfs-to-mount-remote-file-systems-over-ssh) (currently not implemented due to performance and stability issues)  

### Homepage
* Select from available transcriptomes   
* Provide basic information on Trinity and Trinotate versions, analysis pipeline, etc.    

### Data presentation
#### Metadata
* Experimental design  
* Analysis methods (Rmarkdown)  

#### Transcript information  
* Select from available transcriptomes  
* Allow filtering based on keyword/transcript/gene  
* Download selected sequences or annotation table  

#### Expression
* Show selected transcript expression (in a table and graph)  
* Allow consolidating multiple transcripts under one gene name/symbol   

## CFX-qpcR
A shiny server to analyse qPCR results from Biorad CFX platform.  

### Main features
* Select genes/targets  
* Warn about outliers/extreme variability among replicates  
* Display calibration curves plots  
* Save linear regression, R<sup>2</sup> and calculated *Eff* in a table available for download  
