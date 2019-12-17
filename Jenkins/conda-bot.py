import re
import sys
import yaml

from packaging import version
import conda.cli.python_api
from github import Github


TOKEN = sys.argv[1]
g = Github(TOKEN)
repo = g.get_repo('vanheeringen-lab/snakemake-workflows')
open_pulls = list(repo.get_pulls(state='open'))

# ignore packages in the ignore file
with open('./envs/ignore.txt') as f:
    ignore = f.read().splitlines()

for file in repo.get_contents("envs", ref='develop'):
    updated = dependencies = False
    newfile = ''

    pull_request_title = f'auto-update: {file.name}: '

    lines = list(filter(None, file.decoded_content.decode("utf-8").split('\n')))
    yaml_dict = yaml.safe_load(file.decoded_content.decode("utf-8"))

    for line in lines:
        if not any(unwanted in line for unwanted in ['pip', 'github']) and \
               any(dependency in line for dependency in yaml_dict['dependencies'] if isinstance(dependency, str)):

            channel, package, old_version = filter(None, re.split(':| |=', yaml.safe_load(line)[0]))

            # check for newer version
            stdout, stderr, return_code = conda.cli.python_api.run_command('search', package, channel=channel)

            versions = []
            for subline in stdout.split('\n'):
                if not subline == '':
                    output = [sub for sub in subline.split(' ') if sub != '']

                    if output[-1] == channel:
                        versions.append(output)

            if len(versions):
                new_version = versions[-1][1]

                # check if package was updated, and whether or not we want automatic updates for it
                if version.parse(new_version) > version.parse(old_version) and \
                        package not in ignore:

                    # replace the old version with the new version
                    line = line.replace(old_version, new_version)

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
            file = repo.get_file_contents(file.path, ref='develop')
            sb = repo.get_branch(source_branch)
            repo.create_git_ref(ref='refs/heads/' + target_branch, sha=sb.commit.sha)
            repo.update_file(f'envs/{file.name}', f"auto-update: {file.path}", str.encode(newfile), file.sha, branch=target_branch)

            repo.create_pull(title=pull_request_title,head=target_branch, base=source_branch,
                             body='This is an automated pull-request... :robot: ')

            break
