#!/usr/bin/env python
"""
This is the user's entry-point for the seq2science tool.
"""
import sys
import argparse
import os
import shutil
import webbrowser
import re
import inspect
import yaml
import contextlib


# we need to be able to get the parser from the file without a valid seq2science installation
try:
    import snakemake
    import seq2science

    full_import = True
except ImportError:
    full_import = False

try:
    import mamba
    conda_frontend = "mamba"
except ImportError:
    conda_frontend = "conda"



def main():
    # set helpful paths
    base_dir = os.path.dirname(inspect.getfile(seq2science))
    workflows_dir = os.path.join(base_dir, "workflows")

    parser = seq2science_parser(workflows_dir)
    args = parser.parse_args()

    # now run the command
    if args.command == "init":
        dir_path = args.dir if os.path.isabs(args.dir) else os.path.join(os.getcwd(), args.dir)
        _init(args, workflows_dir, dir_path)
    elif args.command == "run":
        config_path = args.configfile if os.path.isabs(args.configfile) else os.path.join(os.getcwd(), args.configfile)
        _run(args, base_dir, workflows_dir, config_path)
    elif args.command == "explain":
        config_path = args.configfile if os.path.isabs(args.configfile) else os.path.join(os.getcwd(), args.configfile)
        _explain(args, base_dir, workflows_dir, config_path)
    elif args.command == "clean":
        _clean(base_dir)
    elif args.command == "docs":
        _docs()


def seq2science_parser(workflows_dir="./seq2science/workflows/"):
    """
    Make the seq2science parser.
    """
    # setup the parser
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-v", "--version", action="version", version=f"seq2science: v{seq2science.__version__}"
    )
    subparsers = parser.add_subparsers(dest="command")
    subparsers.required = True
    init = subparsers.add_parser(
        "init",
        help="Initialise a workflow with an example config and samples file.",
        description="Each workflow requires a configuration and samples file to run. "
                    'Running "seq2science init {workflow}" initialises a default '
                    "configuration and samples file for the specific workflow."
    )
    global run
    run = subparsers.add_parser(
        "run",
        help="Run a complete workflow.",
        description="Run a complete workflow. This requires that a config and samples file "
        "are either present in the current directory, or passed as an argument."
    )
    explain = subparsers.add_parser(
        "explain",
        help="Write a materials & methods section.",
        description="Explains what has/will be done for the workflow. This prints a string which can serve"
                    " as a skeleton for your material & methods section."
    )
    clean = subparsers.add_parser(
        "clean",
        help="Remove all cached sample layouts and conda environments.",
        description="At the start of each workflow run, seq2science starts with installing environments for each "
                    "rule. It also stores the layout of public samples in its cache. These environments can get large "
                    "and it might be best to remove them when you are done with an analysis. \n"
                    "seq2science clean will clean up these files for you."
    )
    docs = subparsers.add_parser(
        "docs",
        description="The docs command tries to open your browser and open the docs' webpage, "
        "if that didn't work it prints the url.",
        help="Take me to the docs!",
    )

    # both init and run can use all workflows
    for subparser in [init, run, explain]:
        subparser.add_argument("workflow", choices=[dir.replace("_", "-") for dir in os.listdir(workflows_dir)])

    # setup init arguments
    init.add_argument(
        "--dir",
        default=".",
        metavar="PATH",
        help="The path to the directory where to initialise the config and samples files.",
    )

    run.add_argument(
        "-c",
        "--configfile",
        default="./config.yaml",
        metavar="FILE",
        help="The path to the config file.",
    )
    global core_arg
    core_arg = run.add_argument(
        "-j", "--cores",
        metavar="N",
        type=int,
        required=True,
        help="Use at most N cores in parallel. Must be at least 2.",
    )
    run.add_argument(
        "-n", "--dryrun",
        help="Do not execute anything, and display what would be done.",
        action='store_true'
    )
    run.add_argument(
        "--unlock",
        help="Remove a lock on the working directory.",
        action='store_true'
    )
    run.add_argument(
        "--snakemakeOptions",
        nargs="+",
        action=_StoreDictKeyPair,
        metavar="KEY=VAL",
        help="Extra arguments to pass along to snakemake. An example could be seq2science run "
        "alignment --cores 12 --snakemakeOptions resources={mem_gb:100} local_cores=3. Here we pass local_cores as "
        "KEY=VALUE and additional resources can even be passed along in a dictionary. Take a look at the snakemake API "
        "for a complete list of all possible options: "
        "https://snakemake.readthedocs.io/en/latest/api_reference/snakemake.html",
    )
    run.add_argument(
        "--profile",
        metavar="PROFILE NAME",
        help="Use a snakemake/seq2science profile. Profiles can be taken from: https://github.com/snakemake-profiles",
    )
    run.add_argument(
        "-k", "--keep-going",
        help="Go on with independent jobs if a job fails.",
        action='store_true'
    )
    run.add_argument(
        "--rerun-incomplete",
        help="Re-run all jobs the output of which is recognized as incomplete.",
        action='store_true'
    )

    explain.add_argument(
        "-c",
        "--configfile",
        default="./config.yaml",
        metavar="FILE",
        help="The path to the config file.",
    )

    return parser


def _init(args, workflows_dir, config_path):
    """
    Initialise a config.yaml and samples.tsv from the relevant workflow.
    """
    for file in ["samples.tsv", "config.yaml"]:
        src = os.path.join(workflows_dir, args.workflow.replace("-", "_"), file)
        dest = os.path.join(os.path.dirname(config_path), file)

        copy_file = True
        if os.path.exists(dest):
            choices = {"yes": True, "y": True, "no": False, "n": False}

            sys.stdout.write(
                f"File: {dest} already exists. Do you want to overwrite it? (yes/no) "
            )
            while True:
                choice = input().lower()
                if choice in choices:
                    copy_file = choices[choice]
                    break
                else:
                    print("Please respond with yes (y) or no (n).")

        if copy_file:
            shutil.copyfile(src, dest)


def _run(args, base_dir, workflows_dir, config_path):
    """
    Run a complete workflow.
    """
    if not os.path.exists(config_path):
        sys.stdout.write(
            f"The config file: {config_path} does not exist.\nProvide a path to the config file with "
            f"--config or if you do not have a config file run:\n"
            f"seq2science init {args.workflow}\n"
        )
        sys.exit(1)

    # parse the args
    parsed_args = {"snakefile": os.path.join(workflows_dir, args.workflow.replace("-", "_"), "Snakefile"),
                   "cores": args.cores,
                   "use_conda": True,
                   "conda_frontend": conda_frontend,
                   "conda_prefix": os.path.join(base_dir, ".snakemake"),
                   "dryrun": args.dryrun,
                   "keepgoing": args.keep_going,
                   "unlock": args.unlock,
                   "force_incomplete": args.rerun_incomplete,
                   "overwrite_threads": core_parser(args.cores)}

    # get the additional snakemake options
    snakemake_options = args.snakemakeOptions if args.snakemakeOptions is not None else dict()
    snakemake_options.setdefault("config", {}).update({"rule_dir": os.path.join(base_dir, "rules")})
    # parse the profile
    snakemake_options["configfiles"] = [config_path]
    if args.profile is not None:
        profile_file = snakemake.get_profile_file(args.profile, "config.yaml")
        if profile_file is None:
            print("Error: profile given but no config.yaml found.")
            sys.exit(1)
        snakemake_options["configfiles"] += [profile_file]
        profile = yaml.safe_load(open(profile_file).read())
        if "cores" in profile and parsed_args["cores"] is None:
            parsed_args["cores"] = profile["cores"]
    if "overwrite_threads" in snakemake_options:
        parsed_args["overwrite_threads"].update(snakemake_options["overwrite_threads"])

    parsed_args.update(snakemake_options)

    parsed_args["cores"] = int(parsed_args["cores"])

    # run snakemake
    exit_code = snakemake.snakemake(**parsed_args)
    sys.exit(0) if exit_code else sys.exit(1)


def _explain(args, base_dir, workflows_dir, config_path):
    """
    Run a complete dryrun workflow, then return the explanations of each rule used.
    """
    if not os.path.exists(config_path):
        sys.stdout.write(
            f"The config file: {config_path} does not exist.\nProvide a path to the config file with "
            f"--config or if you do not have a config file run:\n"
            f"seq2science init {args.workflow}\n"
        )
        sys.exit(1)

    # parse the args
    # snakemake_options.setdefault("config", {}).update({"rule_dir": os.path.join(base_dir, "rules")})
    parsed_args = {"snakefile": os.path.join(workflows_dir, args.workflow.replace("-", "_"), "Snakefile"),
                   "cores": 999,
                   "dryrun": True,
                   "forceall": True,
                   "quiet": False,
                   "config": {"rule_dir": os.path.join(base_dir, "rules"),
                              "explain_rule": True},
                   "configfiles": [config_path]}

    rules_used = {"start": f"\nPreprocessing of reads was done automatically with workflow tool "
                           f"seq2science v{seq2science.__version__} (https://doi.org/10.5281/zenodo.3921913)."}

    def log_handler(log):
        if log["level"] == "job_info" and "msg" in log and log["msg"] is not None and log["name"] not in rules_used:
            rules_used[log["name"]] = log["msg"]

    parsed_args["log_handler"] = [log_handler]

    # run snakemake (silently)
    with open(os.devnull, "w") as null:
        with contextlib.redirect_stdout(null), contextlib.redirect_stderr(null):
            exit_code = snakemake.snakemake(**parsed_args)

    print(" ".join(rules_used.values()))
    sys.exit(0) if exit_code else sys.exit(1)


def _clean(base_dir):
    """
    Clean the .snakemake folder and layouts cache.
    """
    # remove the snakemake cache
    shutil.rmtree(os.path.join(base_dir, ".snakemake"), ignore_errors=True)

    # remove seq2science caches
    shutil.rmtree(os.path.expanduser(f'~/.config/seq2science'), ignore_errors=True)

    print("All cleaned up!")


def _docs():
    """
    Open a webbrowser to the docs, if that fails simply print the url.
    """
    url = "https://vanheeringen-lab.github.io/seq2science"
    if not webbrowser.open(url):
        print(url)


class _StoreDictKeyPair(argparse.Action):
    """
    Class that allows us to take key=value pairs from command-line to feed to
    snakemake. Solution taken from:
    https://stackoverflow.com/questions/29986185/python-argparse-dict-arg
    """

    def __init__(self, option_strings, dest, nargs=None, **kwargs):
        self._nargs = nargs
        super(_StoreDictKeyPair, self).__init__(
            option_strings, dest, nargs=nargs, **kwargs
        )

    def __call__(self, parser, namespace, values, option_string=None):
        # TODO: cleanup
        my_dict = {}
        for kv in values:
            k, v = kv.split("=")

            if ":" in v:
                pair = list(filter(None, re.split('{|:| |}', v)))
                assert len(pair) == 2
                if pair[1].lower() == 'true':
                    pair[1] = True
                v = {pair[0]: int(pair[1]) if isinstance(pair[1], str) and pair[1].isdigit() else pair[1]}
            elif "[" in v:
                v = re.sub("\[|\]", "", v).split(",")
            try:
                my_dict[k] = int(v)
            except:
                if k not in my_dict:
                    my_dict[k] = v
                else:
                    my_dict[k].update(v)

        setattr(namespace, self.dest, my_dict)


def core_parser(cores):
    """
    Alignment uses a pipe, and snakemake does not handle scaling that, so
    we do it for snakemake.
    """
    cores = int(cores)
    sorters = ["samtools_presort"]
    aligners = ["bwa_mem", "bowtie2_align", "hisat2_align", "star_align"]

    d_sorters_threads = 2
    d_aligner_threads = 10
    desired_threads = d_sorters_threads + d_aligner_threads
    if cores <= 1:
        run.print_help()
        try:
            raise argparse.ArgumentError(core_arg, "specify at least two cores.")
        except argparse.ArgumentError as e:
            print()  # empty line
            print(e)
            sys.exit(1)

    elif cores < desired_threads:
        # scale the threads accordingly
        d_sorters_threads = max([1, round(d_sorters_threads / desired_threads * cores)])
        d_aligner_threads = cores - d_sorters_threads

    overwrite_threads = dict()
    for aligner in aligners:
        overwrite_threads[aligner] = d_aligner_threads
    for sorter in sorters:
        overwrite_threads[sorter] = d_sorters_threads

    return overwrite_threads
