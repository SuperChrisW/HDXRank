import argparse
from pepGraph_embedding import generate_embedding


def main(input_data, output_path):
    pass

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Python project main entry.")
    parser.add_argument("--input", type=str, required=True, help="Input data to process.")
    parser.add_argument("--output", type=str, required=True, help="Output file path to save results.")
    args = parser.parse_args()


    ### generate_embedding(root_dir, task_table, protein_name=None, N_model=1)
    ### preprocessing()
    ### prediction()
    main(args.input, args.output)