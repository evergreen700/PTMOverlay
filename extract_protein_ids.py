from pyteomics import pepxml
import glob
import os
import multiprocessing

file_dir = "file/path" #replace file/path with the path to the folder of .pepXML files
acetyl_files = glob.glob(os.path.join(file_dir, "**", "*.pepXML"), recursive=True)

protein_ids = set()

def process_file(file):
    """Extract unique proteins from a .pepXML file (excluding 'rev' proteins)."""
    proteins = set()
    try:
        for entry in pepxml.read(file):
            protein = entry['search_hit'][0]['proteins'][0]['protein']
            if not protein.startswith("rev"):
                proteins.add(protein)
    except Exception as e:
        print(f"Error processing {file}: {e}")
    return proteins

if __name__ == "__main__":
    ctx = multiprocessing.get_context("forkserver") 
    with ctx.Pool(processes=multiprocessing.cpu_count()) as pool:
        results = pool.map(process_file, acetyl_files)

    # Merge results
    for result in results:
        protein_ids.update(result)

    print(protein_ids)
