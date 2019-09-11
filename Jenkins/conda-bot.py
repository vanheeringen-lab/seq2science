import re
import sys

from packaging import version
import conda.cli.python_api
from github import Github


TOKEN = sys.argv[1]
g = Github(TOKEN)
repo = g.get_repo('vanheeringen-lab/snakemake-workflows')
open_pulls = list(repo.get_pulls(state='open'))


for file in repo.get_contents("envs", ref='develop'):
    updated = dependencies = False
    newfile = ''

    pull_request_title = f'auto-update: {file.name}: '

    lines = list(filter(None, file.decoded_content.decode("utf-8").split('\n')))
    for line in lines:
        # skip lines until the dependencies are listed
        if not dependencies:
            if 'dependencies:' in line:
                dependencies = True
        else:
            # ignore lines starting with pip, and deeper indented lines
            if 'pip:' not in line and not line.startswith('    -'):
                # parse the line
                channel, package, old_version = list(filter(None, list(filter(None, re.split(':| |=|\n', ''.join(re.split('-', line, maxsplit=1)))))))

                # check for newer version
                stdout, stderr, return_code = conda.cli.python_api.run_command('search', package, channel=channel)
                most_recent = stdout.split('\n')[-2]  # last line is empty line, so -2
                new_version = most_recent.split()[1]

                # check if package was updated, and whether or not we want automatic updates for it
                if version.parse(new_version) > version.parse(old_version) and \
                        package not in ['python']:

                    # replace the old version with the new version
                    line = f'  - {channel}::{package}={new_version}'

                    # add to title
                    pull_request_title += f'{package}={old_version} -> {new_version}; '
                    updated = True

        newfile += f"{line}\n"

    pull_request_title = pull_request_title[:-2]

    if updated:
        source_branch = 'develop'
        target_branch = f'auto_{file.name}'

        still_open = any([f'auto-update: {file}:' in pull_request.title
                          for pull_request in open_pulls]) or \
                     target_branch in [branch.name for branch in repo.get_branches()]

        if not still_open:
            file = repo.get_file_contents(file.path)
            sb = repo.get_branch(source_branch)
            repo.create_git_ref(ref='refs/heads/' + target_branch, sha=sb.commit.sha)
            repo.update_file(f'envs/{file.name}', f"auto-update: {file}", str.encode(newfile), file.sha, branch=target_branch)

            repo.create_pull(title=pull_request_title,head=target_branch, base=source_branch,
                             body='This is an automated pull-request... :robot: ')


