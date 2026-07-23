import os
import re
import glob

def update_rust_file(filepath):
    with open(filepath, 'r') as f:
        content = f.read()

    # Regex patterns for Rust constants
    n_stars_pattern = re.compile(r'const\s+N_STARS\s*:\s*usize\s*=\s*\d+\s*;')
    dt_pattern = re.compile(r'const\s+DT\s*:\s*f32\s*=\s*[0-9.eE-]+f32\s*;|const\s+DT\s*:\s*f32\s*=\s*[0-9.eE-]+\s*;')
    steps_pattern = re.compile(r'const\s+STEPS\s*:\s*usize\s*=\s*\d+\s*;')

    new_content = content
    
    if n_stars_pattern.search(new_content):
        new_content = n_stars_pattern.sub('const N_STARS: usize = 5000;', new_content)
        print(f"  Updated N_STARS in {filepath}")
        
    if dt_pattern.search(new_content):
        new_content = dt_pattern.sub('const DT: f32 = 0.0001;', new_content)
        print(f"  Updated DT in {filepath}")
        
    if steps_pattern.search(new_content):
        new_content = steps_pattern.sub('const STEPS: usize = 1000000;', new_content)
        print(f"  Updated STEPS in {filepath}")

    if new_content != content:
        with open(filepath, 'w') as f:
            f.write(new_content)
        return True
    return False

def update_python_file(filepath):
    with open(filepath, 'r') as f:
        content = f.read()

    # Regex patterns for Python variables
    n_stars_pattern = re.compile(r'^\s*N_STARS\s*=\s*\d+', re.MULTILINE)
    dt_pattern = re.compile(r'^\s*DT\s*=\s*[0-9.eE-]+', re.MULTILINE)
    steps_pattern = re.compile(r'^\s*STEPS\s*=\s*\d+', re.MULTILINE)

    new_content = content

    if n_stars_pattern.search(new_content):
        new_content = n_stars_pattern.sub('N_STARS = 5000', new_content)
        print(f"  Updated N_STARS in {filepath}")

    if dt_pattern.search(new_content):
        new_content = dt_pattern.sub('DT = 0.0001', new_content)
        print(f"  Updated DT in {filepath}")

    if steps_pattern.search(new_content):
        new_content = steps_pattern.sub('STEPS = 1000000', new_content)
        print(f"  Updated STEPS in {filepath}")

    if new_content != content:
        with open(filepath, 'w') as f:
            f.write(new_content)
        return True
    return False

def main():
    nbody_dir = '/home/stephane/Science/GAIA/stitch/Nbody'
    
    print("Updating Rust simulation files...")
    rust_files = glob.glob(os.path.join(nbody_dir, '**/src/*.rs'), recursive=True)
    for filepath in rust_files:
        update_rust_file(filepath)

    print("\nUpdating Python simulation files...")
    python_files = glob.glob(os.path.join(nbody_dir, '*.py'))
    for filepath in python_files:
        update_python_file(filepath)

    print("\nFinished updates!")

if __name__ == '__main__':
    main()
