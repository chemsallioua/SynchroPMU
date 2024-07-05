import sys

def process_header(file_path, num_chanls):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    with open(file_path, 'w') as file:
        for line in lines:
            if line.strip().startswith("// #define NUM_CHANLS"):
                line = f"#define NUM_CHANLS {num_chanls}\n"
                print(f"Updated line: {line.strip()}")
            file.write(line)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: process_header.py <file_path> <num_chanls>")
        sys.exit(1)

    file_path = sys.argv[1]
    num_chanls = sys.argv[2]

    print(f"Processing header file: {file_path} with NUM_CHANLS = {num_chanls}")
    process_header(file_path, num_chanls)
    print("Header file processed successfully.")
