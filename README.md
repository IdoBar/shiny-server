# shinotate
A shiny server app to deliver processed and annotated transcriptome data stored in [Trinotate](https://trinotate.github.io/) db

## Main features
* Accessible from a web browser (from anywhere)
* Provide basic summarized information (sequences, annotation, expression)
* Provide basic plotting (expression) 

## Implementation
### Hosting
* Shiny server hosted on a cloud-service (AWS, NECTAR, DigitalOcean, etc.) - see [details](http://deanattali.com/2015/05/09/setup-rstudio-shiny-server-digital-ocean/)
* Mount USC's main servers using `sshfs` - see [details](https://www.digitalocean.com/community/tutorials/how-to-use-sshfs-to-mount-remote-file-systems-over-ssh) 

### Homepage
* Select from available transcriptomes 
* Provide basic information

### Data presentation
#### Metadata
* Experimental design
* Analysis methods (Rmarkdown)
#### Transcript information
* Allow filtering based on keyword/transcript/gene  
* Download selected sequences or annotation table

#### Expression
* Show selected transcript expression (in a table and a graph)
* Allow consolidating multiple transcripts under one gene name/symbol 
