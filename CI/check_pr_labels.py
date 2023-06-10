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

def whatch_list():
    whatchlist_files = dict()
    whatchlist_files['CI/physmon/reference/performance_ambi_orthogonal.root'] = ['Changes Performance - Ambiguity resolution']
    whatchlist_files['CI/physmon/reference/performance_ambi_seeded.root'] = ['Changes Performance - Ambiguity resolution']
    whatchlist_files['CI/physmon/reference/performance_amvf_orthogonal_hist.root'] = ['Changes Performance - Vertex']
    whatchlist_files['CI/physmon/reference/performance_amvf_seeded_hist.root'] = ['Changes Performance - Vertex']
    whatchlist_files['CI/physmon/reference/performance_amvf_truth_estimated_hist.root'] = ['Changes Performance - Vertex']
    whatchlist_files['CI/physmon/reference/performance_amvf_truth_smeared_hist.root'] = ['Changes Performance - Vertex']
    whatchlist_files['CI/physmon/reference/performance_ckf_orthogonal.root'] = ['Changes Performance - Finding']
    whatchlist_files['CI/physmon/reference/performance_ckf_seeded.root'] = ['Changes Performance - Finding']
    whatchlist_files['CI/physmon/reference/performance_ckf_truth_estimated.root'] = ['Changes Performance - Finding']
    whatchlist_files['CI/physmon/reference/performance_ckf_truth_smeared.root'] = ['Changes Performance - Finding']
    whatchlist_files['CI/physmon/reference/performance_gsf.root'] = ['Changes Performance - Fitting']
    whatchlist_files['CI/physmon/reference/performance_ivf_orthogonal_hist.root'] = ['Changes Performance - Vertex']
    whatchlist_files['CI/physmon/reference/performance_ivf_seeded_hist.root'] = ['Changes Performance - Vertex']
    whatchlist_files['CI/physmon/reference/performance_ivf_truth_estimated_hist.root'] = ['Changes Performance - Vertex']
    whatchlist_files['CI/physmon/reference/performance_ivf_truth_smeared_hist.root'] = ['Changes Performance - Vertex']
    whatchlist_files['CI/physmon/reference/performance_seeding_orthogonal.root'] = ['Changes Performance - Seeding']
    whatchlist_files['CI/physmon/reference/performance_seeding_seeded.root'] = ['Changes Performance - Seeding']
    whatchlist_files['CI/physmon/reference/performance_seeding_truth_estimated.root'] = ['Changes Performance - Seeding']
    whatchlist_files['CI/physmon/reference/performance_truth_tracking.root'] = ['Changes Performance - Fitting']
    whatchlist_files['CI/physmon/reference/tracksummary_ckf_orthogonal_hist.root'] = ['Changes Performance - Finding']
    whatchlist_files['CI/physmon/reference/tracksummary_ckf_seeded_hist.root'] = ['Changes Performance - Finding']
    whatchlist_files['CI/physmon/reference/tracksummary_ckf_truth_estimated_hist.root'] = ['Changes Performance - Finding']
    whatchlist_files['CI/physmon/reference/tracksummary_ckf_truth_smeared_hist.root'] = ['Changes Performance - Finding']
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
    pull.create_issue_comment("this is a test of a comment")

    files = pull.get_files()
    print('List of modified files:')
    for el in files:
        print(f"   * {el.filename}")

    required_labels = set()
    for el in files:
        current_required_labels = whatchlist_files.get(el.filename, [])
        if len(current_required_labels) == 0:
            continue
        for label in current_required_labels:
            required_labels.add(label)

    missing_labels = set()
    for required_label in required_labels:
        if required_label in list_labels:
            continue
        missing_labels.add(required_label)

    if len(missing_labels) != 0:
        print(f"Available Labels: {list_labels}")
        print(f"Required Labels: {required_labels}")
        print(f"Missing Labels: {missing_labels}")
        pull.create_issue_comment(f"This PR modifies the performance of a component but it does not contain the required label. Please add: {missing_labels}")
        raise Exception(f"This PR modifies the performance of a component but it does not contain the required label")

if __name__ == "__main__":
    main()
