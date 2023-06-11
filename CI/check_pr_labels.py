from github import Github
import argparse
import re
import os

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--pull_id', required=True, type=int, nargs=1,
                        help='ID of the PR')
    parser.add_argument('--repository', required=True, type=str, nargs=1,
                        help='name of the repository, e.g. "acts-project/acts"')
    return parser.parse_args()

def whatch_list():
    whatchlist_files = dict()
    whatchlist_files['^.github/.*'] = ['Infrastructure']
    whatchlist_files['^CI/.*'] = ['Infrastructure']
    whatchlist_files['^cmake/.*'] = ['Infrastructure']
    whatchlist_files['^Core/.*'] = ['Component - Core']
    whatchlist_files['^Examples/.*'] = ['Component - Examples']
    whatchlist_files['^Fatras/.*'] = ['Component - Fatras']
    whatchlist_files['^Plugins/.*'] = ['Component - Plugins']
    whatchlist_files['^docs'] = ['Component - Documentation']
#    whatchlist_files['^CI/physmon/reference.*'] = ['Changes Performance']
#    whatchlist_files['^CI/physmon/reference/.*_ambi_.*.root'] = ['Changes Performance - Ambiguity resolution']
#    whatchlist_files['^CI/physmon/reference/.*_amvf_.*.root'] = ['Changes Performance - Vertex']
#    whatchlist_files['^CI/physmon/reference/.*_ivf_.*.root'] = ['Changes Performance - Vertex']
#    whatchlist_files['^CI/physmon/reference/.*_ckf_.*.root'] = ['Changes Performance - Finding']
#    whatchlist_files['^CI/physmon/reference//*_seeding_/*.root'] = ['Changes Performance - Seeding']
#    whatchlist_files['^CI/physmon/reference/performance_gsf.root'] = ['Changes Performance - Fitting']
#    whatchlist_files['^CI/physmon/reference/performance_truth_tracking.root'] = ['Changes Performance - Fitting']
    return whatchlist_files

def main():       
    args = parse_arguments()
    github_repository_name = args.repository if type(args.repository) == str else args.repository[0]
    pull_id = args.pull_id if type(args.pull_id) == int else args.pull_id[0]

    print(f'Checking labels for PR #{pull_id} from project: {github_repository_name}')
    
    list_labels = set()
    whatchlist_files = whatch_list()

    github_token = ""
    try:
        github_token = str(os.getenv('GITHUB_TOKEN'))
    except Exception:
        raise Exception("Env variables are not properly set! Check the .env file is present and/or the env variables are set.")
    
    github = Github(github_token)
    repository = github.get_repo(github_repository_name)

    # from https://pygithub.readthedocs.io/en/latest/github_objects/PullRequest.html
    pull = repository.get_pull(pull_id)

    labels = pull.get_labels()
    for label in labels:
        list_labels.add(label.name)

    print('This PR is marked with the following labels:', list_labels)
    #pull.create_issue_comment("this is a test of a comment")
    #pull.add_to_labels("enhancement")

    files = pull.get_files()
    print('List of modified files:')
    for el in files:
        print(f" * {el.filename}")

    print('Checking labels ...')
    for pattern in whatchlist_files:
        for el in files:
            if not re.search(pattern, el.filename):
                continue
            # get the required labels for this pattern
            toadd_labels = whatchlist_files[pattern]
            # add label to PR if missing
            for label in toadd_labels:
                if label not in list_labels:
                    print(f" * adding label: {label}")
                    list_labels.add(label)
                    outcome = pull.add_to_labels(label)
                    if outcome != 200:
                        raise Exception(f'Could not add label "{lavel}" to PR')
            break

    if 'Changes Performance' in list_labels:
        pull.create_issue_comment(f"This PR modifies the performance of a component")
        #raise Exception(f"This PR modifies the performance of a component but it does not contain the required label")

if __name__ == "__main__":
    main()
