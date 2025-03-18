import os
import re
import mmap
from concurrent.futures import ProcessPoolExecutor

# Compile regex patterns once
RMSD_PATTERN = re.compile(r'RMSD: ([+-]?\d*\.\d+(?:[eE][+-]?\d+)?)')
ERROR_PATTERN = re.compile(r'maximum measured error: ([+-]?\d*\.\d+(?:[eE][+-]?\d+)?)')

def extract_data(filename):
    """ Efficiently reads a file and extracts RMSD, error, and coefficients """
    try:
        with open(filename, 'r') as f, mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ) as mm:
            content = mm.read().decode()

        # Extract values using regex
        coefficients = []
        lines = content.splitlines()
        for i, line in enumerate(lines):
            if i < 6:
                coefficients.append(line)
        
        rmsd_match = RMSD_PATTERN.search(content)
        error_match = ERROR_PATTERN.search(content)

        if rmsd_match and error_match:
            return filename, float(rmsd_match.group(1)), float(error_match.group(1)), coefficients
    except Exception as e:
        return None  # Skip unreadable/corrupt files

def find_best_files(folder_path, num_workers=8):
    """ Uses parallel processing to find the best RMSD and error files """
    best_rmsd, best_error = float('inf'), float('inf')
    best_rmsd_file, best_error_file = None, None
    best_rmsd_coefs, best_error_coefs = None, None

    files = [os.path.join(folder_path, f) for f in os.listdir(folder_path) if os.path.isfile(os.path.join(folder_path, f))]

    # Parallel file processing
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        for result in executor.map(extract_data, files):
            if result:
                filename, rmsd, max_error, coefs = result

                if rmsd < best_rmsd:
                    best_rmsd, best_rmsd_file, best_rmsd_coefs = rmsd, filename, coefs

                if max_error < best_error:
                    best_error, best_error_file, best_error_coefs = max_error, filename, coefs

    return (best_rmsd_file, best_rmsd, best_rmsd_coefs), (best_error_file, best_error, best_error_coefs)

def main(folder_path):
    best_rmsd_result, best_error_result = find_best_files(folder_path)

    if best_rmsd_result[0]:
        print(f"Lowest RMSD: {best_rmsd_result[1]} from file {best_rmsd_result[0]}")
        for coef in best_rmsd_result[2]:
            print(f"{coef}")

    if best_error_result[0]:
        print(f"Lowest maximum measured error: {best_error_result[1]} from file {best_error_result[0]}")
        for coef in best_error_result[2]:
            print(f"{coef}")

if __name__ == "__main__":
    folder_path = "results"
    main(folder_path)
