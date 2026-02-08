from pathlib import Path

def fa_name(syst: str):
  return "output/" + syst + "_fixedalpha.txt"

def conv_name(syst: str):
  return "output/" + syst + "_convergence.txt"

def lists_to_file(filepath, *lists, headers=None):
  """
  Write multiple lists to a file as columns.
  
  Args:
    filepath: Path to output file
    *lists: Variable number of lists to write as columns
    headers: Optional list of column headers
  """
  filepath = Path(filepath)
  
  if not lists:
    return
  
  # Determine number of rows (length of longest list)
  max_length = max(len(lst) for lst in lists)
  
  # Use provided headers or generate default ones
  if headers is None:
    headers = [f"col_{i}" for i in range(len(lists))]
  
  with open(filepath, 'w') as f:
    # Write header
    f.write(' '.join(headers) + '\n')
    
    # Write data rows
    for i in range(max_length):
      row = [str(lst[i]) if i < len(lst) else '' for lst in lists]
      f.write(' '.join(row) + '\n')