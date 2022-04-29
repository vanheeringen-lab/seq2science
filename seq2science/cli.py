#!/usr/bin/env python
"""
This is the user's entry-point for the seq2science tool.
"""
import os
import re
import sys
import argparse
import argcomplete
import shutil
import inspect

# we need to be able to get the parser from the file without a valid seq2science installation
try:
    import seq2science
    from seq2science.logging import log_welcome
except ImportError:
    pass


def _import():
    """
    this function serves that we can do imports as late as possible, for faster auto-completion
    """
    global webbrowser, contextlib, yaml, psutil, snakemake, logger, setup_logger, xdg, datetime, _logging
    import yaml
    import psutil
    import webbrowser
    import contextlib
    import datetime

    import snakemake
    from snakemake.logging import logger, setup_logger, _logging
    import xdg


def seq2science_main():
    # set helpful paths
    base_dir = os.path.dirname(inspect.getfile(seq2science))
    workflows_dir = os.path.join(base_dir, "workflows")

    parser = seq2science_parser(workflows_dir)
    args = parser.parse_args()

    # most imports after argparsing for faster tab-completion
    _import()

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


def deseq2science_main():
    # set helpful paths
    base_dir = os.path.dirname(inspect.getfile(seq2science))

    parser = deseq2science_parser()
    args = parser.parse_args()

    # most imports after argparsing for faster tab-completion
    _import()

    # now run the command
    _deseq(args, base_dir)


def seq2science_parser(workflows_dir="./seq2science/workflows/"):
    """
    Make the seq2science parser.
    """
    # setup the parser
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--version", action="version", version=f"seq2science: v{seq2science.__version__}")
    subparsers = parser.add_subparsers(dest="command")
    subparsers.required = True
    init = subparsers.add_parser(
        "init",
        help="Initialise a workflow with an example config and samples file.",
        description="Each workflow requires a configuration and samples file to run. "
        'Running "seq2science init {workflow}" initialises a default '
        "configuration and samples file for the specific workflow.",
    )
    global run
    run = subparsers.add_parser(
        "run",
        help="Run a complete workflow.",
        description="Run a complete workflow. This requires that a config and samples file "
        "are either present in the current directory, or passed as an argument.",
    )
    explain = subparsers.add_parser(
        "explain",
        help="Write a materials & methods section.",
        description="Explains what has/will be done for the workflow. This prints a string which can serve"
        " as a skeleton for your material & methods section.",
    )
    clean = subparsers.add_parser(
        "clean",
        help="Remove all cached sample files and conda environments.",
        description="At the start of each workflow run, seq2science starts with installing environments for each "
        "rule. It also stores the GEO soft files of public samples in its cache. These environments can get"
        " large and it might be best to remove them when you are done with an analysis. \n"
        "seq2science clean will clean up these files for you.",
    )
    docs = subparsers.add_parser(
        "docs",
        description="The docs command tries to open your browser and open the docs' webpage, "
        "if that didn't work it prints the url.",
        help="Take me to the docs!",
    )

    # init, run and explain can use all workflows
    for subparser in [init, run, explain]:
        subparser.add_argument(
            "workflow", metavar="WORKFLOW", choices=[dir.replace("_", "-") for dir in os.listdir(workflows_dir)]
        )

    # init arguments
    init.add_argument(
        "--dir",
        default=".",
        metavar="PATH",
        help="The path to the directory where to initialise the config and samples files.",
    )

    init.add_argument(
        "-f",
        "--force",
        default=False,
        help="Overwrite existing samples.tsv and config.yaml silently.",
        action="store_true",
    )

    global core_arg
    core_arg = run.add_argument(
        "-j",
        "--cores",
        metavar="N",
        type=int,
        # required=True,  # --dryruns and --profile can overwrite None
        help="Use at most N cores in parallel. Must be at least 2. When "
        "executing on a cluster, this number controls the maximum number"
        "of parallel jobs.",
    )
    run.add_argument(
        "-n", "--dryrun", help="Do not execute anything, and display what would be done.", action="store_true"
    )
    run.add_argument("-r", "--reason", help="Print the reason for each executed rule.", action="store_true")
    run.add_argument(
        "--skip-rerun", help="Skip the check if samples or configuration has been changed.", action="store_true"
    )
    run.add_argument("-k", "--keep-going", help="Go on with independent jobs if a job fails.", action="store_true")
    run.add_argument(
        "--rerun-incomplete", help="Re-run all jobs the output of which is recognized as incomplete.", action="store_true"
    )
    run.add_argument("--unlock", help="Remove a lock on the working directory.", action="store_true")
    explain.add_argument("--hyperref", help="Print urls as html hyperref", action="store_true")
    # run/explain arguments
    for subparser in [run, explain]:
        subparser.add_argument(
            "--snakemakeOptions",
            nargs="+",
            action=_StoreDictKeyPair,
            metavar="KEY=VAL",
            help="Extra arguments to pass along to snakemake. An example could be seq2science run "
            "alignment --cores 12 --snakemakeOptions resources={mem_gb:100} local_cores=3. "
            "Here we pass local_cores as KEY=VALUE and additional resources can even be passed along in a dictionary. "
            "Take a look at the snakemake API  for a complete list of all possible options: "
            "https://snakemake-api.readthedocs.io/en/latest/api_reference/snakemake.html",
        )
        global profile_arg
        profile_arg = subparser.add_argument(
            "-p",
            "--profile",
            metavar="PROFILE NAME",
            help="Use a seq2science profile. Profiles can be taken from: https://github.com/s2s-profiles",
        )
        subparser.add_argument(
            "-c",
            "--configfile",
            default="./config.yaml",
            metavar="FILE",
            help="The path to the config file.",
        )
        subparser.add_argument(
            "--debug",
            action="store_true",
            help="""For developers "only": prints helpful error messages to debug issues.""",
        )

    # enable tab completion
    # exclusion only works on the main parser unfortunately, but it's better than nothing,
    # plus it might be supported later?
    argcomplete.autocomplete(parser, exclude=["-c", "-p", "-k" "-r" "-n", "-j", "-h", "-v"])

    return parser


def deseq2science_parser():
    parser = argparse.ArgumentParser(
        description="DESeq2 wrapper that works with s2s samples.tsv and counts.tsv",
    )
    parser.add_argument(
        "-v", "--version", action="version", version=f"seq2science: v{seq2science.__version__}"
    )
    # DESeq2 wrapper arguments
    parser.add_argument(
        "-d",
        "--design",
        default="",  # must be an empty string for R's arg checking + docs CLI argument
        help="design contrast (e.g. column_knockouts_controls)",
    )
    parser.add_argument(
        "-s",
        "--samples",
        default="",  # must be an empty string for R's arg checking + docs CLI argument
        help="samples.tsv with a column containing the contrast arguments "
             "(sample order consistent with the counts.tsv)",
    )
    parser.add_argument(
        "-c",
        "--counts",
        default="",  # must be an empty string for R's arg checking + docs CLI argument
        help="counts.tsv(.gz) (sample order consistent with the samples.tsv)",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        default="",  # must be an empty string for R's arg checking + docs CLI argument
        help="output directory",
    )
    parser.add_argument(
        "-sc",
        "--single-cell",
        action='store_true',
        default=False,
        help="use if the counts are Single Cell data",
    )
    parser.add_argument(
        "--docs",
        action='store_true',
        help="open de DESeq2 wrapper documentation (with examples!)",
    )
    argcomplete.autocomplete(parser)

    return parser


def _init(args, workflows_dir, config_path):
    """
    Initialise a config.yaml and samples.tsv from the relevant workflow.
    """
    for file in ["samples.tsv", "config.yaml"]:
        src = os.path.join(workflows_dir, args.workflow.replace("-", "_"), file)
        dest = os.path.join(os.path.dirname(config_path), file)

        copy_file = True
        if os.path.exists(dest) and args.force is False:
            choices = {"yes": True, "y": True, "no": False, "n": False}

            sys.stdout.write(f"File: {dest} already exists. Do you want to overwrite it? (yes/no) ")
            while True:
                choice = input().lower()
                if choice in choices:
                    copy_file = choices[choice]
                    break
                else:
                    print("Please respond with yes (y) or no (n).")

        if copy_file:
            shutil.copyfile(src, dest)


def subjectively_prettier_error(arg, message):
    """raise an exception and catch it for a subjectively prettier message"""
    try:
        raise argparse.ArgumentError(arg, message)
    except argparse.ArgumentError as err:
        print(f"\n{err}")
        sys.exit(1)


def add_profile_args(profile_file, parsed_args):
    """read profile and add new arguments to parsed args"""
    profile = yaml.safe_load(open(profile_file).read())

    # don't overwrite CLI arguments
    for k, v in profile.items():
        if k not in parsed_args:
            parsed_args[k] = int(v) if isinstance(v, str) and v.isdigit() else v

            # when reading e.g. resources it gets treated as a list, split it into key value pairs
            if isinstance(parsed_args[k], list) and all("=" in item for item in parsed_args[k]):
                parsed_args[k] = {item.split("=")[0]: float(item.split("=")[1]) for item in parsed_args[k]}

        elif k in parsed_args and isinstance(parsed_args[k], dict):
            for k2, v2 in parsed_args[k].items():
                if k2 not in parsed_args[k]:
                    parsed_args[k][k2] = int(v2) if isinstance(v2, str) and v2.isdigit() else v2


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
    parsed_args = {
        "snakefile": os.path.join(workflows_dir, args.workflow.replace("-", "_"), "Snakefile"),
        "use_conda": True,
        "conda_frontend": "mamba",
        "conda_prefix": os.path.join(base_dir, ".snakemake"),
        "dryrun": args.dryrun,
        "printreason": args.reason,
        "keepgoing": args.keep_going,
        "unlock": args.unlock,
        "force_incomplete": args.rerun_incomplete,
    }

    # get the additional snakemake options
    snakemake_options = args.snakemakeOptions if args.snakemakeOptions is not None else dict()
    snakemake_options.setdefault("config", {}).update({"rule_dir": os.path.join(base_dir, "rules")})
    snakemake_options["configfiles"] = [config_path]
    for key, value in snakemake_options.items():
        if not isinstance(value, str):
           continue
        if value.lower() == "true":
           snakemake_options[key] = True
        if value.lower() == "false":
           snakemake_options[key] = False

    parsed_args.update(snakemake_options)

    # parse the profile
    if args.profile is not None:
        profile_file = snakemake.get_profile_file(args.profile, "config.yaml")
        if profile_file is None:
            subjectively_prettier_error(profile_arg, "profile given but no config.yaml found.")
        add_profile_args(profile_file, parsed_args)

    # cores
    if args.cores:  # command-line interface
        parsed_args["cores"] = args.cores
    elif parsed_args.get("cores"):  # profile
        parsed_args["cores"] = int(parsed_args["cores"])
    elif parsed_args["dryrun"]:
        parsed_args["cores"] = 999
    else:
        parsed_args["cores"] = 0

    if parsed_args["cores"] < 2:
        subjectively_prettier_error(core_arg, "specify at least two cores.")

    # when running on a cluster assume cores == nodes (just like snakemake does)
    if "cluster" in parsed_args and not "nodes" in parsed_args:
        parsed_args["nodes"] = parsed_args["cores"]

    # store how seq2science was called
    parsed_args["config"]["cli_call"] = sys.argv

    # core_parser(parsed_args)  # TODO: do we still need this?
    parsed_args["config"].update({"cores": parsed_args["cores"]})
    resource_parser(parsed_args)

    # run snakemake/seq2science
    #   1. pretty welcome message
    setup_seq2science_logger(parsed_args)
    log_welcome(logger, args.workflow)

    if args.debug:
        # dump the parsed args as readable json
        import json
        logger.info(json.dumps(parsed_args, sort_keys=True, indent=2))

    if not args.skip_rerun or args.unlock:
        #   2. start a dryrun checking which files need to be created, and check if
        #      any params changed, which means we have to remove those files and
        #      continue from there
        logger.info(
            "Checking if seq2science was run already, if something in the configuration was changed, and if so, if "
            "seq2science needs to re-run any jobs."
        )

        with seq2science.util.CaptureStdout() as targets, seq2science.util.CaptureStderr() as errors:
            exit_code = snakemake.snakemake(
                **{
                    **parsed_args,
                    **{
                        "list_params_changes": True,
                        "quiet": False,
                        "log_handler": lambda x: None,  # don't show any of the logs
                        "keep_logger": True,
                    },
                }
            )
        if args.debug:
            nl = "\n"
            logger.info(f"""Targets:\n{nl.join(sorted(targets))}\n\n""")
            logger.info(f"""Errors:\n{nl.join(sorted(errors))}\n\n""")
        if not exit_code:
            sys.exit(1)

        #   3. check which files would need a rerun, and exclude files we do
        #      not want to consider:
        #      - genome files, since provider will change to local
        regex_patterns = [
            "(\/.+){2}.*\.(fa(\.fai|.sizes)?|gaps\.bed)$",  # match genome files
            "(\/.+){2}.*\.annotation\.(bed|gtf)$",  # match gene annotations
        ]
        targets = [target for target in targets if not any(re.match(pattern, target) for pattern in regex_patterns)]

        #   4. if there are any targets left, force to recreate those targets plus the final results (rule seq2science)
        if len(targets):
            targets += ["seq2science"]
            parsed_args["forcerun"] = targets
            parsed_args["targets"] = targets
            parsed_args["forcetargets"] = True
            parsed_args["keep_logger"] = True
        logger.info("Done. Now starting the real run.")

    logger.printreason = parsed_args["printreason"]
    logger.stream_handler.setStream(sys.stdout)
    parsed_args["config"]["no_config_log"] = True

    #   5. start the "real" run where jobs actually get started
    exit_code = snakemake.snakemake(**parsed_args)

    #   6. output exit code 0 for success and 1 for failure
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
    parsed_args = {
        "snakefile": os.path.join(workflows_dir, args.workflow.replace("-", "_"), "Snakefile"),
        "dryrun": True,
        "forceall": True,
        "quiet": False,
    }

    # get the additional snakemake options
    snakemake_options = args.snakemakeOptions if args.snakemakeOptions is not None else dict()
    snakemake_options.setdefault("config", {}).update({"rule_dir": os.path.join(base_dir, "rules")})
    snakemake_options["configfiles"] = [config_path]
    parsed_args.update(snakemake_options)

    # parse the profile
    if args.profile is not None:
        profile_file = snakemake.get_profile_file(args.profile, "config.yaml")
        if profile_file is None:
            subjectively_prettier_error(profile_arg, "profile given but no config.yaml found.")
        add_profile_args(profile_file, parsed_args)

    # cores
    parsed_args["cores"] = 999
    parsed_args["config"]["hyperref"] = args.hyperref

    # starting message
    if args.hyperref:
        rules_used = {
            "start": f"Preprocessing of reads was done automatically by "
            f"<a href=https://doi.org/10.5281/zenodo.3921913>seq2science v{seq2science.__version__}</a> "
            f"using the {args.workflow} workflow."
        }
    else:
        rules_used = {
            "start": f"Preprocessing of reads was done automatically by "
            f"seq2science v{seq2science.__version__} (https://doi.org/10.5281/zenodo.3921913) "
            f"using the {args.workflow} workflow."
        }

    def log_handler(log):
        if (
            log["level"] == "job_info"
            and "msg" in log
            and log["name"] not in rules_used
            and log["msg"] not in (list(rules_used.values()) + [None])
        ):
            rules_used[log["name"]] = log["msg"]

    parsed_args["log_handler"] = [log_handler]
    parsed_args["config"]["explain_rule"] = True

    # run snakemake (silently)
    with open(os.devnull, "w") as null:
        with contextlib.redirect_stdout(null), contextlib.redirect_stderr(null):
            success = snakemake.snakemake(**parsed_args)

    if args.debug:
        print(f"Explain output:\n{rules_used}\n\n")

    if success:
        print(" ".join(rules_used.values()))
        sys.exit(0)
    else:
        print(
            "Oh no! Something went wrong... "
            "Please let us know: https://github.com/vanheeringen-lab/seq2science/issues "
        )
        sys.exit(1)


def _clean(base_dir):
    """
    Clean the .snakemake folder and cache.
    """
    # remove the snakemake cache
    shutil.rmtree(os.path.join(base_dir, ".snakemake"), ignore_errors=True)

    # remove seq2science caches
    shutil.rmtree(os.path.expanduser(os.path.join(xdg.XDG_CACHE_HOME, "seq2science")), ignore_errors=True)

    # remove historic seq2science cache location
    shutil.rmtree(os.path.expanduser(f"~/.config/seq2science/"), ignore_errors=True)

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
        super(_StoreDictKeyPair, self).__init__(option_strings, dest, nargs=nargs, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        # TODO: cleanup
        my_dict = {}
        for kv in values:
            k, v = kv.split("=")

            if ":" in v:
                assert "}" in v if "{" in v else True, (
                    f"\n\nThe dictionary you provided on the command line ('{k}') cannot be parsed:"
                    f"\n'{v}' (TIP: is there a space in there perhaps?)\n"
                )
                pair = list(filter(None, re.split("{|:| |}", v)))
                assert len(pair) == 2, (
                    f"\n\nThe dictionary you provided on the command line ('{k}') cannot be parsed:"
                    f"\n'{v}' contains a broken key-value pair: '{pair}' (TIP: is there a space in there perhaps?)\n"
                )
                if pair[1].lower() == "true":
                    pair[1] = True
                elif pair[1].lower() == "false":
                    pair[1] = False
                elif isinstance(pair[1], str) and pair[1].isdigit():
                    pair[1] = int(pair[1])
                v = {pair[0]: pair[1]}
            elif "[" in v:
                v = re.sub("\[|\]", "", v).split(",")
            else:
                assert k != "config", (
                    f"\n\nThe dictionary you provided on the command line ('{k}') cannot be parsed:"
                    f"\n'{v}' should be a key-value pair (TIP: use '{v}:True' perhaps?)\n"
                )

            if k in my_dict:
                my_dict[k].update(v)
            else:
                my_dict[k] = int(v) if isinstance(v, str) and v.isdigit() else v

        setattr(namespace, self.dest, my_dict)


def core_parser(parsed_args):
    """
    Alignment uses a pipe, and snakemake does not handle scaling that, so
    we do it for snakemake.
    """
    cores = parsed_args["cores"]
    sorters = ["samtools_presort"]
    aligners = ["bowtie2_align", "bwa_mem", "bwa_mem2", "hisat2_align", "minimap2_align", "star_align"]

    d_sorters_threads = 2
    d_aligner_threads = 10
    desired_threads = d_sorters_threads + d_aligner_threads
    overwrite_threads = dict()
    # scale if aligners+sorters use more threads than allowed
    if cores < desired_threads:
        # scale the threads accordingly
        d_sorters_threads = max([1, cores * d_sorters_threads // desired_threads])
        d_aligner_threads = max([1, cores - d_sorters_threads])

        for aligner in aligners:
            overwrite_threads[aligner] = d_aligner_threads
        for sorter in sorters:
            overwrite_threads[sorter] = d_sorters_threads

    # finally use the user-specified threads
    if "overwrite_threads" in parsed_args:
        for k, v in parsed_args["overwrite_threads"].items():
            overwrite_threads[k] = v

    # convert to snakemake's key=value format
    overwrite_threads = [f"{k}={v}" for k, v in overwrite_threads.items()]

    parsed_args["overwrite_threads"] = overwrite_threads


def resource_parser(parsed_args):
    """
    Make sure to properly parse the resources.
    Memory limit changes depending on local execution or cluster
    """
    # set some defaults
    parsed_args["resources"] = {
        **{"parallel_downloads": 3, "genomepy_downloads": 1, "deeptools_limit": 16, "R_scripts": 1},
        **parsed_args.get("resources", {}),
    }

    if "mem_mb" in parsed_args["resources"]:
        # convert memory to gigabytes
        parsed_args["resources"]["mem_gb"] = round(parsed_args["resources"]["mem_mb"] / 1024.0)
        del parsed_args["resources"]["mem_mb"]

    # no need to get system limit when specified
    if "mem_gb" in parsed_args["resources"]:
        return

    if "cluster" in parsed_args:
        # if running on a cluster assume no limit on memory (unless specified)
        parsed_args["resources"]["mem_gb"] = 999999
    else:
        # otherwise assume system memory
        mem = psutil.virtual_memory().total / 1024 ** 3
        parsed_args["resources"]["mem_gb"] = round(mem)


def setup_seq2science_logger(parsed_args):
    setup_logger()
    if not parsed_args["dryrun"]:
        seq2science_logfile = os.path.abspath(
            "seq2science." + datetime.datetime.now().isoformat().replace(":", "") + ".log"
        )
        logger.logfile = seq2science_logfile
        logger.get_logfile = lambda: seq2science_logfile
        logger.logfile_handler = _logging.FileHandler(seq2science_logfile)
        logger.logger.addHandler(logger.logfile_handler)


def _deseq(args, base_dir):
    # docs
    if args.docs is True:
        url = "https://vanheeringen-lab.github.io/seq2science/content/DESeq2.html"
        if not webbrowser.open(url):
            print(url)
        return

    # insufficient args
    if not all(len(a) > 0 for a in [args.counts, args.design, args.outdir, args.samples]):
        logger.info("4 arguments expected: contrast, samples_file, counts_file and outdir.")
        return

    import hashlib
    import subprocess as sp

    def conda_path(yml):
        """
        Find the path to a snakemake conda environment.
        Does not work with containers.
        """
        # mimic snakemake's conda env hashing:
        # https://github.com/snakemake/snakemake/blob/
        # 1349254de57ced0466bfe1e98f0f7c2a1b247d2d/snakemake/deployment/conda.py#L91
        env_dir = os.path.dirname(os.path.dirname(yamlfile))
        env_dir = os.path.join(env_dir, ".snakemake")
        env_dir = os.path.realpath(env_dir)
        md5hash = hashlib.md5()
        md5hash.update(env_dir.encode())
        md5hash.update(open(yml, "rb").read())
        env_hash = md5hash.hexdigest()[:8]
        path = os.path.join(env_dir, env_hash)
        return path

    def subprocess_run(_cmd):
        retcode = sp.call(_cmd, shell=True)
        print("")  # no newline otherwise
        if retcode != 0:
            sys.exit(retcode)

    # find/create the deseq2 conda env
    yamlfile = os.path.join(base_dir, "envs", "deseq2.yaml")
    env_prefix = conda_path(yamlfile)
    if not os.path.exists(env_prefix):
        logger.info(f"Creating conda environment seq2science/envs/deseq2.yaml...")
        cmd = f"mamba env create -p {env_prefix} -f {yamlfile} -q > /dev/null"
        subprocess_run(cmd)

    # we don't even need to activate the env
    rscript = os.path.join(env_prefix, "bin", "Rscript")
    script = os.path.join(base_dir, "scripts", "deseq2", "deseq2.R")
    cmd = f"{rscript} {script} {args.design} {args.samples} {args.counts} {args.outdir} {args.single_cell}"
    subprocess_run(cmd)

    # example command:
    #
    # deseq2science \
    # -d batch+condition_day2_day0 \
    # -s tests/deseq2/rna/samples.tsv \
    # -c tests/deseq2/rna/counts/GRCh38.p13-counts.tsv \
    # -o tests/local_test_results/deseq2science
