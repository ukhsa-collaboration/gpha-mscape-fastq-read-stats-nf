# gpha-mscape-fastq-read-stats-nf

## What is this?

A nextflow pipeline wrapping [gpha-mscape-fastq-read-stats](https://github.com/ukhsa-collaboration/gpha-mscape-fastq-read-stats),
which is an argparse tool for generating basic statistics, given a fastq.gz file as input, designed to run on [mSCAPE](https://mscape.climb.ac.uk/).

## How do I use this?

On the CLIMB infrastructure, you'd run a command not dissimilar to the following, replacing `<CLIMB-ID>` with an actual CLIMB ID.

```bash
nextflow run                      \
	main.nf                       \
	--profile docker              \
	--unique_id <CLIMB-ID>        \
	-e.ONYX_DOMAIN=$ONYX_DOMAIN   \
    -e.ONYX_TOKEN=$ONYX_TOKEN
```

The site nextflow config stored at `/etc/nextflow.config` is required as a `-c` argument if your `nextflow` command has not already been aliased to include this.

## What's the theory behind it?

Something like this:

```mermaid
graph TD

graphStart(["Start"])
userInput{"Parse commandline args"}
onyxQuery["Retrieve record from Onyx save as samplesheet"]
endiannessDecision{"Single or paired end?"}
preprocessingOption1["Retrieve and stage unclassified.fastq.gz"]
preprocessingOption2["Retreive and stage unclassified.1.fastq.gz and unclassified.2.fastq.gz"]
processing["Process files sequentially"]
output["Result published to fastq_read_stats_run/<CLIMB-ID>"]
graphEnd(["End"])

graphStart --> userInput
userInput -- CLIMB ID<br>provided via args ---> onyxQuery 
userInput -- "Samplesheet<br>provided via args" ---> endiannessDecision 
onyxQuery --> endiannessDecision
endiannessDecision -- Single ---> preprocessingOption1
endiannessDecision -- Paired ---> preprocessingOption2
preprocessingOption1 & preprocessingOption2 --> processing
processing --> output --> graphEnd
```
