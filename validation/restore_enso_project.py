"""Stage recovered ENSO inputs. Does not run R or rebuild climate inputs."""
from pathlib import Path
import gzip, shutil

root=Path(__file__).resolve().parent
project=root/'enso_project'
targets={'data_qc': [('inputs/yield_qc.csv','soybean_Brazil_Argentina_analysis_ready_yield.csv')],
         'enso_indices':[('inputs/enso_by_year.csv','ENSO_indices_by_soybean_harvest_year_1974_2024.csv')],
         'metadata':[('inputs/spatial_units_calendars.csv','soybean_spatial_units_with_calendars.csv')]}
for folder,pairs in targets.items():
    dest=project/folder;dest.mkdir(parents=True,exist_ok=True)
    for source,name in pairs:shutil.copy2(root/source,dest/name)
merged=project/'data_merged';merged.mkdir(parents=True,exist_ok=True)
with gzip.open(root/'inputs/recovered_enso_local_panel.csv.gz','rb') as source:
    with (merged/'soybean_ENSO_spatial_panel_1974_2024.csv').open('wb') as target:
        shutil.copyfileobj(source,target)
scripts=project/'scripts';scripts.mkdir(exist_ok=True)
for source in (root/'archived_R_scripts').glob('*.R'):shutil.copy2(source,scripts/source.name)
print('Recovered ENSO project staged at:',project)
print('The local climate panel and matching seven-index outputs are not staged.')
