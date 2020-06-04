# Snakemake help 

> This document contains useful instructions for bioinformaticians to better understand snakemake and how we make use of it in the Huntsman Lab.

#### Asli Munzur, April 2nd, 2020

This README is intended to help colleagues who haven't used snakemake before. We walk you through how it is designed, and give a few eaxmples of how we can various types of files through it.  

You can refer <a href="https://snakemake.readthedocs.io/en/stable/" target="_blank">here</a> to read about it, but I will summarize below how we adopted it for our needs. 

The idea is to be able to automate certain repetitive parts of the analysis. For example, regardless of the sample we analyze, we will need to make an sce from the data (or perform qc on it). Being able to automate it saves us time when we have many samples to analyze.  

We divide the snakefile into 'rules'. These are the types of analysis we want to run. For example, making an sce would be a rule, performing qc would be another rule. 

A rule looks like this: 

```
rule some_rule: 
  input: 
    path_to_input_file 
  output: 
    path_to_output_file
  params: 
    any special parameters you need in your script
  shell: 
    path_to_script
    declare parameters 
```

For example, here is the rule we use to make an sce from the raw counts: 

```
rule make_sce:
  input:
    "../data/raw/{id}"
  output:
    "../data/processed/{id}/sce.rds"
  params:
    "{id}"
  shell:
    "Rscript pipeline/data_preparation/convert_to_sces.R \
     --path_to_10X {input} \
     --output_file_name {output} \
     --id {params} "
```

Did {id} catch your eye? That is a wildcard. Wildcards are like variables, they can take various values. That is one of the powers of snakemake. I can define a list of sample ids to be supplied to the wildcards shown in the rule above, and the rule will be repeated for each of the wildcards, meaning the script you intend to run will be run for each of the ids supplied. Let's take a look at a slightly more complicated rule: 

```
rule batch_correction:
  params:
    curr_dir = os.getcwd(), 
    ids_integration = ids_integration, 
    id_type = id_type
  input:
    sce_clus
  output:
    report = batch_correction_report,
    output_uncorrected = with_batch,
    output_corrected = without_batch
  shell:
    "Rscript -e \"rmarkdown::render('pipeline/batch_normalization/batch_normalization.Rmd',\
     output_file='{params.curr_dir}/{output.report}', \
     knit_root_dir='{params.curr_dir}',\
     params = list(ids ='{params.ids_integration}', id_type = '{params.id_type}', output_uncorrected = '{output.output_uncorrected}', output_corrected = '{output.output_corrected}'))\" "
```

We use this rule above to remove batch effect in more than one sample. We have a document named 'Orchestrating Single-Cell Analysis with Bioconductor' in the onboarding page, it has information about batch effects. This rule above is a little more complicated because it is a RMarkdown document, and they usually contain more information in the rule. Now we will go through the components of the rule one by one. This next bit is taken from the yaml header in the markdown document itself, you'll see how they relate to one another:    

```
params:
  ids: 'ids_integration'
  output_uncorrected: 'output_uncorrected'
  output_corrected: 'output_corrected'
  id_type: 'id_type'
```

Basically, this file has some variables like the samples we use and we declare them in the header. When we look at this header, we see that the file uses four variables: ids, output_uncorrected, output_corrected and id_type. We also see this:  

```
ids: 'ids_integration'
```

The variable 'ids' is equal to 'ids_integration'. What this means is that Snakemake passes the string 'ids_integration' for ids. Now, we should go back to the Snakefile and figure out what ids_integration is equal to. This is what we see at the top of the script:  

```
ids_integration = ['DH4-DH17-DH10-DH3-DH15-DH16']
```

What all this means is that our script has a certain unknown, which is filled in by what Snakemake passes to it while running it. This is true for all the variables declared in the rule. The values Snakemake passes for them are given in the Snakefile:  

```
id_type = ['ENDOMETRIAL']
```

However, not all the variables declared in the yaml header of the markdown document has to be declared elsewhere in the Snakefile. Sometimes they are given somewhere in the rule. We see that the variables used in the markdown file named output_uncorrected and output_corrected are equal to without_batch and with_batch, respectively. without_batch and with_batch are defined in the output portion of the rule. Both are file paths, and essentially they are declared at the top of the Snakefile: 

```
# BATCH CORRECTION 
with_batch = expand('../data/batch_corrected/{ids_integration}/uncorrected.rds', ids_integration = ids_integration)
without_batch = expand('../data/batch_corrected/{ids_integration}/corrected.rds', ids_integration = ids_integration) 
```

This is basically the case:  
unknown in the YAML header -> snakemake passes a value for that -> passed value is taken from output -> output takes the value from elsewhere in the Snakefile

#### Now we will look at the 'params' section of the rule:  
**curr_dir = os.getcwd()**  
This just tells the document where to run it. We always add this to all the rules about markdown documents.  
**ids_integration = ids_integration** and **id_type = id_type**  
This tells us that the rule uses a variable called ids_integration, and it is equal to another variable named 'ids_integration' that is defined in the Snakefile. I always make sure that they have the same name, so that I dont have to deal with various names all over the place, but it is a matter of preference. You can give them different names for sure.  

#### Next let's take a look at the output section of the rule:  
This section contains two types of outputs. One is knitted markdown file, and the other is the any kind of outputs we may have in the markdown file. Sometimes our markdown files will also compute outputs (like R objects) that we may want to save for future use. We need to tell R where to save those objects, and this is where the output comes into play. In the rule above, we see that we have three outputs: 

```
output:
 report = batch_correction_report,
 output_uncorrected = with_batch,
 output_corrected = without_batch
```

'report' is where we the knitted html (or pdf, up to you but we go with html) file will be saved. We see that it is equal to batch_correction_report, so we go up the Snakefile to see what batch_correction_report is: 

```
batch_correction_report = expand('../reports/batch_normalization/{ids_integration}/batch_correction.html', ids_integration = ids_integration)
```

Same thing goes for the two other outputs given in the rule. Since batch_correction_report is where the knitted file will be saved, we don't need to call it anywhere in the script. However, the other two outputs are for saving an object for future use, so we do call them in the script when saving, like this:  

```
saveRDS(object, file = output_uncorrected)
```


#### Now we look at the final section, shell:  
Rules that execute markdown files usually have four lines in this part. The first line starting with Rscript is telling the terminal to execute an R script. The rest is just where the files is. You'll notice that the lines end with \, this is just to tell R that the code continues in the next line. When you replicate it, make sure to keep all the punctuation the same. In the next line, you see 'output_file'. This is how we tell R where to download knitted output file. Again, make sure you use the same punctuation! The quotation marks we have here helps R interpret this as a string, which is indeed a string since it is a file path. Don't forget the comma: 

```
output_file='{params.curr_dir}/{output.report}', \
```

Next line in the rule is this: 

```
knit_root_dir='{params.curr_dir}',\
```

Add this to ALL the rules concerning RMarkdown files. It helps R set the directory correctly.  

Next we see this in the rule: 

```
params = list(ids ='{params.ids_integration}', id_type = '{params.id_type}', output_uncorrected = '{output.output_uncorrected}', output_corrected = '{output.output_corrected}'))\" "
```

Now, this a little bit confusing. We use this to pass arguments to the markdown file while we run it. If you have more than one argument to pass, they should be given as a list. Here, we define all the variables we have in the yaml header in the RMarkdown file we wish to execute through Snakemake. More imprtantly, we also tell R where to find them. Note that we didn't just do this: 
```
params = list(ids ='{ids_integration}', id_type = '{id_type}', output_uncorrected = '{output_uncorrected}', output_corrected = '{output_corrected}'))\" "
```

This would be wrong. R needs to know where each of these variables come from because they can come from more than one place. For example, ids_integration and id_type come from 'params' while output_uncorrected and output_corrected come from 'output'. Therefore it is important to distinguish them. 

Congratulations on making it this far! I found that running markdown files through Snakemake is more difficult than running other types of scripts because of the extra setting we have to specify. Now we will go through the other type of script you may have: just regular R scripts. Let's take a look at the beginning of the make_sce script:  

```
parser <- ArgumentParser(description = "Convert 10X output to sce")

parser$add_argument('--path_to_10X', metavar='DIRECTORY', type='character',
                    help="Path to cellranger output directory")

parser$add_argument('--id', metavar='VARIABLE', type='character',
                    help="unique id of the input dataset")

parser$add_argument('--output_file_name', metavar='FILE', type='character',
                    help="Path to raw sce")

args <- parser$parse_args()
```

Before you panic, read a little bit about argument parser. You can follow this pattern in the future scripts you write. First, we define the variables we have in our script. In the case I showed above, we need three variables to run this script:  
- path to the 10X data where we have the gene counts saved  
- the sample id, such as DH22, DH25 ...  
- the path where we want to save the computed sce object  

Basically Snakemake uses argument parser to pass these input to our script. We start by giving a description of what our script is doing. In this case, we want to convert 10X output to an sce. Therefore, we write:  

```
parser <- ArgumentParser(description = "Convert 10X output to sce")
```

Next, we need to add arguments. Follow the same punctuation and add the variables you wish to run in your script. Reading a little bit about argument parser will help you specify options, like metavar and help. When you are done, just do this: 

```
args <- parser$parse_args()
```

This helps R keep the arguments in one place, called R. When you are done and ready to run your script, do this at the end: 

```
convert_to_sces(path_to_10X = args$path_to_10X, 
                output_file_name = args$output_file_name,
                id = args$id)
```


By doing this, you tell R where to find the arguments. Script needs path_to_10X as an input. We just saved it to 'args', so we call it by doing 

```
path_to_10X = args$path_to_10X
```

#### Rule all  
Rule all is another part of Snakemake. There are various ways of using it, but I use to tell Snakemake which objects I want to make.  

```
rule all:
  input:
    sce_raw,
    sce_qc,
    sce_norm,
    sce_clus,
    sce_uncorrected, # only combined sces, without removing batch effects
    seurat_integ, # combined sces with batch effects removed
    integration_report, # dim reduction plots before and after batch effect removal
    separate_clustering_report, # cluster an sce alone
    with_batch,
    without_batch,
    batch_correction_report,
    summary_stats_report,
    edgeR_report,
    paired_dge_basic_report
```

You'll notice that I write my objects and file paths that the rules will need at the top of the Snakefile. Then I add them all to here. Let's say I need to generate the edgeR_report for some samples. I will look at that specific rule to see which wildcards it needs to run. I will fill those wildcards appropriately with the samples I am interested in running, then I will comment out objects I don't need to compute. This is more of a habit honestly. If sce_raw is already computed, Snakemake won't make the object anyways. This is how rule all would look like in this case: 

```
rule all:
  input:
    # sce_raw,
    # sce_qc,
    # sce_norm,
    # sce_clus,
    # sce_uncorrected, # only combined sces, without removing batch effects
    # seurat_integ, # combined sces with batch effects removed
    # integration_report, # dim reduction plots before and after batch effect removal
    # separate_clustering_report, # cluster an sce alone
    # with_batch,
    # without_batch,
    # batch_correction_report,
    # summary_stats_report,
    edgeR_report,
    # paired_dge_basic_report
```

That's it! If you have these basics, you will be able to learn more complicated rules easily. The snakefile has many more rules that you can take a look at, but they all follow these basic rules.  

#### Now we have a few notes about how to run this specific Snakefile:  
It's super easy to run these steps:  
- make sce  
- perform qc  
- normalize  
- dim reduction  
- clustering  

They all use the same wildcards. So lets say you need to run DH4, DH5, DH6 and DH7 through the steps mentioned above. You don't need to run them one by one, you can just fill the id list with the samples ids you are interested in, and all the samples will go through the pipeline one by one. So cool! However, a few of the rules aren't as straightforward. For example, let's say you want to make the integration report, then you need to integrate one by one. Snakefile is changing and getting better everyday, so this possibly will be fixed at some point in the future. Moreover, because the Snakefile may change later on, I pushed a copy of it to Huntsman Lab's Github account where the things I mentioned above are valid. The commit id is ca2be66c53f6423e00f4f9db503a9f0d55c80197. You can pull from there and use the file.  

#### Running Snakemake  
Go to the terminal in the RStudio. Don't use your own terminal! Type "snakemake -n" (without quotation marks). This is a dry run, it will tell you if you made any syntax errors etc. To actually run snakemake, all you have to do is type "snakemake" to the terminal.  

Happy learning, and ask you colleagues if you have any questions!  

