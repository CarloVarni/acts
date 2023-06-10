from github import Github
import argparse
import os

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--pull_id', required=True, type=int, nargs=1,
                        help='ID of the PR')
    parser.add_argument('--repository', required=True, type=str, nargs=1,
                        help='name of the repository, e.g. "acts-project/acts"')
    return parser.parse_args()
    
def main():       
    args = parse_arguments()
    github_repository_name = args.repository if type(args.repository) == str else args.repository[0]
    pull_id = args.pull_id if type(args.pull_id) == int else args.pull_id[0]

    print(f'Checking labels for PR #{pull_id} from project: {github_repository_name}')
    
    list_labels = set()

    whatlist_files = dict()
    whatlist_files['CI/physmon/reference/performance_amvf_orthogonal_hist.root'] = ['Changes Performance - Vertex']
    whatlist_files['CI/physmon/reference/performance_amvf_seeded_hist.root'] = []
    whatlist_files['CI/physmon/reference/performance_amvf_truth_estimated_hist.root'] = []
    whatlist_files['CI/physmon/reference/performance_amvf_truth_smeared_hist.root'] = []

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

    files = pull.get_files()
    print('List of modified files:')
    for el in files:
        print(f"   * {el.filename}")
        
    for el in files:
        required_labels = whatlist_files.get(el.filename, [])
        if len(required_labels) == 0:
            continue
        for required_label in required_labels:
            assert (required_label in list_labels), f"This PR modifies the performance of a component but it does not contain the required label: '{required_label}'"

if __name__ == "__main__":
    main()
