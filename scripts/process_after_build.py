import sys

def process_header(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    with open(file_path, 'w') as file:
        for line in lines:
            if line.strip().startswith("#define NUM_CHANLS"):
                print(f"Removed line: {line.strip()}")
                line = "// #define NUM_CHANLS 1\n"
                print(f"Added comment line: {line.strip()}")
            file.write(line)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: process_header.py <file_path>")
        sys.exit(1)

    file_path = sys.argv[1]

    print(f"Processing header file: {file_path} to remove NUM_CHANLS and add comment")
    process_header(file_path)
    print("Header file processed successfully.")
