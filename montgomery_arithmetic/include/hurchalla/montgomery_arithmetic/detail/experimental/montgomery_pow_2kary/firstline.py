import sys

def main():
    if len(sys.argv) != 2:
        print("Usage: python3 script.py <filename>")
        sys.exit(1)

    filename = sys.argv[1]

    try:
        with open(filename, 'r') as file:
            lines = file.readlines()
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
        sys.exit(1)

    # Find start and end markers
    try:
        start_index = next(i for i, line in enumerate(lines) if "OVERALL BEST:" in line)
        end_index = next(i for i, line in enumerate(lines) if "Timings By Test Type:" in line)
    except StopIteration:
        print("Error: Could not find required markers in the file.")
        sys.exit(1)

    # Process lines between the markers
    for line in lines[start_index + 1:end_index]:
        parts = line.strip().split()
        if len(parts) != 7:
            continue  # skip malformed lines
        try:
            third_field = int(parts[2])
        except ValueError:
            continue  # skip lines where the third field isn’t an integer
        if third_field < 6:
            print(line.strip())
            return

    print("No line found where the third field is less than 6.")

if __name__ == "__main__":
    main()
