#!/usr/bin/env python3

# steps required
# 1. read in csv - generate summary stats
# 2. create summary stats dictionary
# 3. create analysis dictionary 
# 4. write out onyx json 

#!/usr/bin/env python3

"""Script and associated functions to take unclassified read QC
results and create an onyx analysis table format json file.
"""

# Imports
import argparse
import logging
import os
import sys
from ast import literal_eval
from pathlib import Path

import pandas as pd
import yaml
from onyx import OnyxConfig, OnyxEnv
from onyx_analysis_helper import onyx_analysis_helper_functions as oa

# Set up onyx config
CONFIG = OnyxConfig(
    domain=os.environ[OnyxEnv.DOMAIN],
    token=os.environ[OnyxEnv.TOKEN],
)


# Functions
# Args and logging
def get_args():
    """Get command line arguments"""
    parser = argparse.ArgumentParser(
        prog="Create read QC analysis table",
        description="""Simple script to generate onyx analysis fields from unclassified read QC results
        """,
    )
    parser.add_argument("--climbid", "-i", type=str, required=True, help="Sample ID")
    parser.add_argument(
        "--server",
        "-s",
        type=str,
        required=True,
        choices=["mscape", "synthscape"],
        help="Specify server code is being run on",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        required=True,
        help="Folder to save downloaded files to",
    )
    parser.add_argument(
        "--pipeline_info",
        "-p",
        type=str,
        required=True,
        help="Comma separated str in format: 'pipeline name, version, homepage'",
    )
    parser.add_argument(
        "--version",
        "-v",
        type=str,
        required=True,
        help="Unclassified read QC version",
    )
    parser.add_argument(
        "--readQC_result",
        "-r",
        type=str,
        required=True,
        help="Argument stating whether unclassified read QC ran successfully",
    )
    parser.add_argument(
        "--readQC_file",
        "-f",
        type=str,
        required=True,
        help="Filepath for read QC results",
    )
    parser.add_argument(
        "--context",
        "-c",
        dest="orange_box_version",
        type=str,
        required=False,
        help="Orange box version used to determine if should run.",
    )
    # Add group options to specify onyx behaviour
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--store-onyx",
        required=False,
        action="store_true",
        help="Use this option to store results as an onyx analysis object for later upload",
    )
    group.add_argument(
        "--test-onyx",
        required=False,
        action="store_true",
        help="Use this option to do a test upload and check for errors before attempting an upload to onyx",
    )
    group.add_argument(
        "--prod-onyx",
        required=False,
        action="store_true",
        help="Use this option to upload results to onyx",
    )

    return parser.parse_args()


def set_up_logger(stdout_file):
    """Creates logger for component - all logging messages go to stdout
    log file, error messages also go to stderr log. If component runs
    correctly, stderr is empty.
    """
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter("[%(asctime)s] %(levelname)s: %(message)s")

    out_handler = logging.FileHandler(stdout_file, mode="a")
    out_handler.setFormatter(formatter)
    logger.addHandler(out_handler)

    return logger


def get_unclassified_results(unclass_data_file: Path) -> dict:
    """Parses the result_data file from unclassified read QC to generate summary statistics.
    Arguments:
        all_data_file -- Path to file containing all unclassified read QC results
    Returns:
        Dictionary of unclassified read QC results
    """

    readqc = pd.read_csv(unclass_data_file, sep = "\t")

    result_dict = readqc.agg(
    {
        "mean_phred_score": ["min", "max", "mean","std","median"],
        "read_length": ["min", "max", "mean","std","median"],
        "gc_content": ["min", "max", "mean","std","median"],
        "compression_ratio": ["min", "max", "mean","std","median"],
        "paired_base_complexity": ["min", "max", "mean","std","median"],
        
    }).to_dict()
    result_dict['unclassified_read_count']=readqc.shape[0]
    result_dict

    return result_dict

def make_pipeline_dict(pipeline_str) -> dict:
    """Create dict of pipeline info from input pipeline info string.
    Arguments:
        pipeline_str -- Comma separated string of pipeline info in format
                        'pipeline_name,pipeline_version,pipeline_homepage'
    Returns:
        pipeline_dict -- Dictionary containing pipeline information
    """

    pipeline_list = pipeline_str.split(",")
    pipeline_dict = {
        "name": pipeline_list[0],
        "version": pipeline_list[1],
        "homePage": pipeline_list[2],
    }

    return pipeline_dict

## What else to change here
def create_analysis_fields(
    record_id: str,
    unclassified_version: str,
    unclassified_qc_results: dict,
    server: str,
    pipeline_info: dict,
    orange_box_version: str,
) -> tuple[oa.OnyxAnalysis, int]:
    onyx_analysis = oa.OnyxAnalysis()
    onyx_analysis.add_analysis_details(
        analysis_name="gpha-mscape-fastq-read-stats-nf",
        analysis_description="This is an analysis to perform read QC for unclassified reads",
    )
    onyx_analysis.pipeline_name = pipeline_info["name"]
    onyx_analysis.pipeline_version = pipeline_info["version"]
    onyx_analysis.pipeline_url = pipeline_info["homePage"]
    onyx_analysis.outputs = "No s3 outputs"  # ty:ignore[invalid-assignment]
    # Add settings and onyx analysis hash
    methods_fail = onyx_analysis.add_methods(
        methods_dict={
            "unclassified_version": unclassified_version,
        }
    )

    _, onyx_versions, query_exitcode = oa.get_data_and_versions_from_onyx(
        sample_id=record_id, server=server, fields=["climb_id"]
    )

    if query_exitcode == 0:
        # Add onyx versions and orange box version
        methods_versions_fail = onyx_analysis.add_versions_to_methods(
            tool_versions={"orange_box_version": orange_box_version},
            onyx_versions=onyx_versions,  # onyx versions hash automatically generated
        )
    else:
        methods_versions_fail = False

    headline_result = unclassified_qc_results["unclassified_read_count"]
    results_fail = onyx_analysis.add_results(
        top_result=headline_result, results_dict=unclassified_qc_results
    )
    onyx_analysis.add_server_records(sample_id=record_id, server_name=server)
    required_field_fail, attribute_fail = onyx_analysis.check_analysis_object(
        publish_analysis=False
    )

    if any(
        [
            methods_fail,
            methods_versions_fail,
            results_fail,
            required_field_fail,
            attribute_fail,
        ]
    ):  # noqa SIM108
        exitcode = 1
    else:
        exitcode = 0

    return onyx_analysis, exitcode


# Main script
def main():
    "Main function to process unclassified read QC results ready for submission to onyx"
    args = get_args()
    exitcode = 0

    # Set up log file
    log_file = Path(args.output) / f"{args.climbid}.unclassified.readQC.analysis_fields.log.txt"
    set_up_logger(log_file)

    orange_box_version = (
        args.orange_box_version if args.orange_box_version else "unknown"
    )

    if args.readQC_result == "True":

        # Get information from unclassified read QC result files
        result_dict = get_unclassified_results(args.readQC_file)

    else:

        result_dict = {"unclassified_read_count": "No data - empty file"}
    # Convert pipeline info to dict
    pipeline_dict = make_pipeline_dict(args.pipeline_info)

    # Create onyx analysis dict
    onyx_analysis, exitcode = create_analysis_fields(
        record_id=args.climbid,
        unclassified_version=args.version,
        unclassified_qc_results=result_dict,
        server=args.server,
        pipeline_info=pipeline_dict,
        orange_box_version=orange_box_version,
    )

    # Exit if analysis object not made correctly
    if exitcode == 1:
        logging.error(
            "Invalid attribute in analysis fields submitted, check logs for details"
        )
        return exitcode

    # Store analysis table info in json file
    if args.store_onyx:
        onyx_json_file = (
            Path(args.output) / f"{args.climbid}.unclassified_readQC.analysis_fields.json"
        )
        result_file = onyx_analysis.write_analysis_to_json(result_file=onyx_json_file)
        logging.info("Onyx analysis fields written to file %s", result_file)
        exitcode = 0
        return exitcode

    # Test push data to onyx
    if args.test_onyx:
        result, exitcode = onyx_analysis.write_analysis_to_onyx(
            server=args.server, dryrun=True, publish_analysis=False
        )

    # Push analysis table to onyx
    if args.prod_onyx:
        result, exitcode = onyx_analysis.write_analysis_to_onyx(
            server=args.server, dryrun=False, publish_analysis=False
        )

    return exitcode


if __name__ == "__main__":
    sys.exit(main())
