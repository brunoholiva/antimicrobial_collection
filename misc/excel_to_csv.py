import pandas as pd

def excel_to_csv(excel_file, sheet_name, csv_file, skiprows=0):
    df = pd.read_excel(excel_file, sheet_name=sheet_name, skiprows=skiprows)
    df.to_csv(csv_file, index=False)


excel_to_csv(excel_file='data/raw/stokes_2020.xlsx', 
             sheet_name='S1B', 
             csv_file='data/raw/stokes_2020.csv', 
             skiprows=1)

excel_to_csv(excel_file='data/raw/swanson_2024_data.xlsx', 
             sheet_name=0, 
             csv_file='data/raw/swanson_2024.csv',
             skiprows=1)

excel_to_csv(excel_file='data/raw/wong_2024.xlsx',
             sheet_name='S. aureus growth inhibition',
             csv_file='data/raw/wong_2024.csv')