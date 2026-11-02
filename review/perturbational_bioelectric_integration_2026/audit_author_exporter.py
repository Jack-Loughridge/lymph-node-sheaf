"""Actual downloaded author source-code schema audit; not physiological payload QC."""
from pathlib import Path
import re,json,hashlib
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
def literal_exports(source):
 active='\n'.join(line for line in source.splitlines() if not line.lstrip().startswith('%'))
 return sorted(set(re.findall(r'ret\.(\w+)\s*=',active)))
def main():
 base=ROOT/'data/perturbational_bioelectric_integration_2026';p=base/'astro_python_conversion_aqua_extract.m';s=p.read_text();fields=literal_exports(s)
 result={'evidence_type':'ACTUAL_AUTHOR_CODE_NOT_PHYSIOLOGICAL_ARRAYS','source_url':'https://github.com/poskanzerlab/astro-network-encoding/blob/v1.0/python/conversion/aqua_extract.m','pinned_git_blob':'8ceb1189bc62bf1d637cb76fbc08529297157e52','cached_UTF8_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'literal_export_fields':fields,'event_marks_dynamically_copied':True,'uncaging_frame_presets_copied_from_export_options':True,'raw_laser_voltage_log_exported':False,'continuous_fluorescence_trace_exported':False,'all_cell_catalog_or_masks_exported':False,'zero_event_cells_safely_distinguishable_from_unobserved':False,'cell_identity_is_local_event_membership_not_global_cell_type':True,'exported_centroids_are_event_footprints_not_independent_cell_ROI_geometry':True,'entire12GB_MAT_schema_verified':False,'physiological_files_downloaded':0,'DATA_GO':False}
 (HERE/'AUTHOR_EXPORT_SCHEMA_AUDIT.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
 sources=[{'source_file':p.name,'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'kind':'author_code_not_biological_observations'} for p in sorted(base.glob('astro_*'))]
 (HERE/'RETRIEVED_CODE_INVENTORY.json').write_text(json.dumps(sources,indent=2)+'\n')
if __name__=='__main__':main()
