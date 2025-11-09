# remove_lines.py

import sys

def remove_lines(input_filename, search_string, output_filename):
    """
    Reads lines from input_filename and writes to output_filename
    all lines that do not contain search_string.
    """
    with open(input_filename, 'r', encoding='utf-8') as infile, \
         open(output_filename, 'w', encoding='utf-8') as outfile:
        
        for line in infile:
            if search_string not in line:
                outfile.write(line)

def main():
    # Expect exactly three command-line arguments
    if len(sys.argv) != 4:
        print("Usage: python remove_lines.py <input_file> <search_string> <output_file>")
        sys.exit(1)

    input_filename = sys.argv[1]
    search_string = sys.argv[2]
    output_filename = sys.argv[3]

    remove_lines(input_filename, search_string, output_filename)
    print(f"Lines containing '{search_string}' have been written to '{output_filename}'.")

if __name__ == "__main__":
    main()