# AFM AGRFORMET-D-26-03447 — execução completa da revisão
# Pode executar o arquivo inteiro por Source ou colar todo o código no console.

project_dir <- Sys.getenv("ENSO_SOYBEAN_PROJECT_DIR", unset = "enso_project")
revision_dir <- Sys.getenv("ENSO_REVISION_DIR", unset = getwd())

if (!dir.exists(project_dir)) stop("Pasta do projeto não encontrada.")
if (!dir.exists(revision_dir)) stop("Pasta da revisão não encontrada.")
project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)
revision_dir <- normalizePath(revision_dir, winslash = "/", mustWork = TRUE)
climate_dir <- file.path(project_dir, "climate_mechanisms")

necessarios <- c(
  file.path(revision_dir, c(
    "verify_revision.py", "revision_models.py",
    "run_water_memory.py", "run_climate_revision.py"
  )),
  file.path(revision_dir, "inputs", c(
    "enso_by_year.csv", "published_field_oof.csv"
  )),
  file.path(climate_dir, c(
    "NASA_POWER_grid_manifest.csv",
    "soybean_spatial_units_with_calendars.csv",
    "soybean_climate_indices_grid_year.csv.gz",
    "soybean_water_balance_primary_AWC100.csv"
  ))
)
ausentes <- necessarios[!file.exists(necessarios)]
if (length(ausentes)) {
  stop("Arquivos ausentes:\n", paste(ausentes, collapse = "\n"))
}
for (pasta in c("NASA_POWER_daily_grid", "NASA_POWER_daily_ET0_inputs")) {
  if (!dir.exists(file.path(climate_dir, pasta))) {
    stop("Pasta meteorológica ausente: ", pasta)
  }
}

# A auditoria de TODAS as grades ocorre no bloco Python abaixo.
# Estima somente falhas isoladas de radiação fora de TODOS os ciclos de soja.
# Cada estimativa preserva backup e registro; outras lacunas são listadas juntas.

# Configura Python somente se ainda não estiver iniciado nesta sessão.
if (!requireNamespace("reticulate", quietly = TRUE)) {
  install.packages("reticulate", repos = "https://cloud.r-project.org")
}
if (utils::packageVersion("reticulate") < "1.41.0") {
  stop("Atualize reticulate e reinicie o R antes de executar este arquivo.")
}
if (!reticulate::py_available(initialize = FALSE)) {
  Sys.setenv(
    RETICULATE_PYTHON = "managed",
    OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1"
  )
  reticulate::py_require(
    c("numpy==2.3.5", "pandas==2.2.3", "matplotlib>=3.8,<4"),
    python_version = ">=3.11,<3.14"
  )
}

# Define os caminhos diretamente no mesmo bloco Python que executa os testes.
codigo_python <- paste(
  sprintf("project_dir = %s", encodeString(project_dir, quote = '"')),
  sprintf("revision_dir = %s", encodeString(revision_dir, quote = '"')),
  "\"\"\"Audit every required POWER cache before the water-memory calculation.",
  "",
  "Only isolated, internally bounded radiation gaps outside every assigned crop",
  "cycle can be interpolated. Missingness must not exceed 0.1% of radiation days",
  "in a grid. Backups, before/after hashes and estimates are retained. No other",
  "weather variable is imputed, and no meteorology is downloaded here.",
  "\"\"\"",
  "from pathlib import Path",
  "from datetime import datetime, timezone",
  "import hashlib, json, os, shutil, tempfile",
  "import numpy as np",
  "import pandas as pd",
  "",
  "",
  "def preflight_daily_weather(project_dir, output_dir):",
  "    from run_water_memory import midpoint",
  "",
  "    project_dir, output_dir = Path(project_dir), Path(output_dir)",
  "    output_dir.mkdir(parents=True, exist_ok=True)",
  "    climate_dir = project_dir/'climate_mechanisms'",
  "    manifest = pd.read_csv(climate_dir/'NASA_POWER_grid_manifest.csv',",
  "                           dtype={'climate_location_id': str, 'grid_id': str})",
  "    calendars = pd.read_csv(climate_dir/'soybean_spatial_units_with_calendars.csv',",
  "                            dtype={'climate_location_id': str})",
  "    month_cols = ['sowing_start_month', 'sowing_end_month',",
  "                  'harvest_start_month', 'harvest_end_month']",
  "    mapping = manifest[['climate_location_id', 'grid_id']].drop_duplicates()",
  "    units = calendars[['climate_location_id', 'calendar_group']+month_cols].merge(",
  "        mapping, on='climate_location_id', validate='one_to_one', how='left')",
  "    if units.grid_id.isna().any():",
  "        raise ValueError('Some crop calendars have no grid mapping; inspect metadata.')",
  "    groups = units[['grid_id', 'calendar_group']+month_cols].drop_duplicates()",
  "    if groups.duplicated(['grid_id', 'calendar_group']).any():",
  "        raise ValueError('Ambiguous crop-calendar definition for a POWER grid.')",
  "    expected = pd.date_range('1984-01-01', '2024-07-31', freq='D')",
  "    specifications = [",
  "        ('meteorology', 'NASA_POWER_daily_grid', '.csv.gz',",
  "         ['T2M', 'T2M_MAX', 'T2M_MIN', 'T2MDEW', 'WS2M', 'PRECTOTCORR']),",
  "        ('ET0_inputs', 'NASA_POWER_daily_ET0_inputs', '_ET0.csv.gz',",
  "         ['ALLSKY_SFC_SW_DWN', 'PS'])",
  "    ]",
  "    gaps, files, repairs = [], [], []",
  "    grids = sorted(groups.grid_id.unique())",
  "    stamp = datetime.now(timezone.utc).isoformat()",
  "    for number, grid in enumerate(grids, 1):",
  "        if number == 1 or number % 20 == 0 or number == len(grids):",
  "            print(f'Auditoria meteorologica {number}/{len(grids)}: {grid}', flush=True)",
  "        active = np.zeros(len(expected), dtype=bool)",
  "        for g in groups[groups.grid_id == grid].itertuples(index=False):",
  "            for year in range(1985, 2025):",
  "                start = midpoint(year, g.sowing_start_month, g.sowing_end_month)",
  "                end = midpoint(year, g.harvest_start_month, g.harvest_end_month)",
  "                if end <= start:",
  "                    raise ValueError(f'Invalid crop dates: {grid}/{g.calendar_group}/{year}')",
  "                active |= (expected >= start) & (expected <= end)",
  "        for family, folder, suffix, variables in specifications:",
  "            source = climate_dir/folder/(grid+suffix)",
  "            base = dict(grid_id=grid, family=family, source_file=str(source))",
  "            try:",
  "                table = pd.read_csv(source, float_precision='round_trip')",
  "                missing_cols = set(['date']+variables)-set(table.columns)",
  "                if missing_cols:",
  "                    raise ValueError('Missing columns: '+', '.join(sorted(missing_cols)))",
  "                dates = pd.to_datetime(table.date, errors='coerce')",
  "                if dates.isna().any() or dates.duplicated().any():",
  "                    raise ValueError('Invalid or duplicate daily dates')",
  "                positions = pd.Series(np.arange(len(table)), index=pd.DatetimeIndex(dates))",
  "                aligned = table[variables].copy()",
  "                aligned.index = pd.DatetimeIndex(dates)",
  "                aligned = aligned.reindex(expected)",
  "                missing_count = 0",
  "                file_repairs = []",
  "                for variable in variables:",
  "                    values = pd.to_numeric(aligned[variable], errors='coerce').to_numpy(dtype=float)",
  "                    missing = ~np.isfinite(values) | (values <= -999)",
  "                    missing_count += int(missing.sum())",
  "                    sparse_radiation = variable == 'ALLSKY_SFC_SW_DWN' and missing.mean() <= .001",
  "                    for j in np.flatnonzero(missing):",
  "                        day = expected[j]",
  "                        reason = 'Other meteorological variable; no automatic imputation'",
  "                        can_repair = False",
  "                        if variable == 'ALLSKY_SFC_SW_DWN':",
  "                            if not sparse_radiation:",
  "                                reason = 'Radiation missingness exceeds 0.1% of expected days'",
  "                            elif active[j]:",
  "                                reason = 'Within at least one assigned soybean crop cycle'",
  "                            elif j == 0 or j == len(expected)-1:",
  "                                reason = 'Boundary gap; no two-sided interpolation'",
  "                            elif day not in positions.index:",
  "                                reason = 'Daily row missing; pressure may also be unavailable'",
  "                            elif missing[j-1] or missing[j+1]:",
  "                                reason = 'Consecutive radiation gaps; no automatic interpolation'",
  "                            elif values[j-1] < 0 or values[j+1] < 0:",
  "                                reason = 'Invalid neighboring radiation values'",
  "                            else:",
  "                                can_repair = True",
  "                                reason = 'Single-day linear interpolation outside all crop cycles'",
  "                        record = dict(**base, variable=variable, date=str(day.date()),",
  "                                      within_crop=bool(active[j]), status='unresolved', reason=reason,",
  "                                      old_value=None if not np.isfinite(values[j]) else float(values[j]),",
  "                                      estimated_value=None)",
  "                        if can_repair:",
  "                            record['estimated_value'] = float((values[j-1]+values[j+1])/2)",
  "                            record['neighbor_date_before'] = str(expected[j-1].date())",
  "                            record['neighbor_date_after'] = str(expected[j+1].date())",
  "                            record['neighbor_value_before'] = float(values[j-1])",
  "                            record['neighbor_value_after'] = float(values[j+1])",
  "                            file_repairs.append((int(positions.loc[day]), record))",
  "                        gaps.append(record)",
  "                before = hashlib.sha256(source.read_bytes()).hexdigest()",
  "                after = before",
  "                if file_repairs:",
  "                    changed = table.copy()",
  "                    for row, record in file_repairs:",
  "                        changed.loc[row, record['variable']] = record['estimated_value']",
  "                    backup = Path(str(source)+'.before_radiation_QC.bak')",
  "                    if not backup.exists():",
  "                        with source.open('rb') as src, backup.open('xb') as dst:",
  "                            shutil.copyfileobj(src, dst)",
  "                    descriptor, name = tempfile.mkstemp(prefix='POWER_QC_', suffix='.csv.gz', dir=source.parent)",
  "                    os.close(descriptor)",
  "                    temporary = Path(name)",
  "                    try:",
  "                        changed.to_csv(temporary, index=False, compression='gzip')",
  "                        check = pd.read_csv(temporary, float_precision='round_trip')",
  "                        pd.testing.assert_frame_equal(changed, check, check_dtype=False,",
  "                                                      check_exact=False, rtol=1e-12, atol=1e-12)",
  "                        os.replace(temporary, source)",
  "                    finally:",
  "                        if temporary.exists():",
  "                            temporary.unlink()",
  "                    after = hashlib.sha256(source.read_bytes()).hexdigest()",
  "                    for row, record in file_repairs:",
  "                        record.update(status='interpolated', backup_file=str(backup),",
  "                                      sha256_before=before, sha256_after=after,",
  "                                      timestamp_utc=stamp, unit='MJ/m2/day')",
  "                        repairs.append(record.copy())",
  "                    history = output_dir/'radiation_repair_history.csv'",
  "                    pd.DataFrame([r.copy() for _, r in file_repairs]).to_csv(",
  "                        history, mode='a', header=not history.exists(), index=False)",
  "                    print(f'  {grid}: {len(file_repairs)} estimativa(s) de radiacao registrada(s).', flush=True)",
  "                files.append(dict(**base, missing_before=missing_count,",
  "                                  repaired=len(file_repairs), unresolved=missing_count-len(file_repairs),",
  "                                  sha256_before=before, sha256_after=after, status='audited'))",
  "            except Exception as exc:",
  "                gaps.append(dict(**base, variable='FILE', date='', within_crop=None,",
  "                                 status='unresolved', reason=f'{type(exc).__name__}: {exc}'))",
  "                files.append(dict(**base, status='file_error', error=str(exc)))",
  "    columns = ['grid_id', 'family', 'source_file', 'variable', 'date',",
  "               'within_crop', 'status', 'reason', 'old_value', 'estimated_value']",
  "    gap_table = pd.DataFrame(gaps)",
  "    if gap_table.empty:",
  "        gap_table = pd.DataFrame(columns=columns)",
  "    unresolved = gap_table[gap_table.status == 'unresolved'].copy()",
  "    gap_table.to_csv(output_dir/'all_missing_cells.csv', index=False)",
  "    unresolved.to_csv(output_dir/'lacunas_nao_resolvidas.csv', index=False)",
  "    pd.DataFrame(files).to_csv(output_dir/'file_audit.csv', index=False)",
  "    summary = dict(timestamp_utc=stamp, grid_count=len(grids), file_count=len(files),",
  "                   new_radiation_estimates=len(repairs), unresolved_records=len(unresolved),",
  "                   policy='Only isolated one-day radiation gaps outside ALL assigned crop cycles; missingness <=0.1%; no other variable imputed',",
  "                   previous_repair_logs=[str(p) for p in (output_dir.parent/'data_repairs').glob('*.csv')],",
  "                   passed=unresolved.empty)",
  "    (output_dir/'preflight_summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')",
  "    if not unresolved.empty:",
  "        raise ValueError(f'Auditoria completa: {len(unresolved)} pendencia(s) em '",
  "                         f'{unresolved.grid_id.nunique()} grade(s). Consulte: '",
  "                         f'{output_dir / \"lacunas_nao_resolvidas.csv\"}')",
  "    print(f'Auditoria concluida: {len(grids)} grades; '",
  "          f'{len(repairs)} estimativa(s) nova(s); nenhuma lacuna pendente.', flush=True)",
  "    return summary",
  "",
  "",
  "from pathlib import Path",
  "import sys, runpy, zipfile, traceback, os",
  "import numpy, pandas, matplotlib",
  "",
  "root = Path(revision_dir)",
  "sys.path.insert(0, str(root))",
  "jobs = [",
  "    ('verify_revision.py', []),",
  "    ('run_water_memory.py', ['--project', project_dir, '--output', str(root/'results_water')]),",
  "    ('run_climate_revision.py', ['--project', project_dir, '--water', str(root/'results_water'/'water_initialization_local.csv.gz'), '--output', str(root/'results_climate')])",
  "]",
  "log = root/'EXECUCAO_STATUS.txt'",
  "temporary_zip = root/'RESULTADOS_PARA_DEVOLVER.tmp.zip'",
  "original_argv = list(sys.argv)",
  "log.write_text('RUNNING\\n', encoding='utf-8')",
  "",
  "try:",
  "    print('Python:', sys.version.split()[0], flush=True)",
  "    print('NumPy:', numpy.__version__, '| pandas:', pandas.__version__, flush=True)",
  "    print('Iniciando auditoria completa dos arquivos meteorologicos.', flush=True)",
  "    preflight_daily_weather(Path(project_dir), root/'results_water'/'preflight')",
  "    with log.open('a', encoding='utf-8') as stream:",
  "        stream.write('DONE: meteorological preflight\\n')",
  "    for script, arguments in jobs:",
  "        print('\\nExecutando:', script, flush=True)",
  "        with log.open('a', encoding='utf-8') as stream:",
  "            stream.write('START: '+script+'\\n')",
  "        sys.argv = [script] + arguments",
  "        runpy.run_path(str(root/script), run_name='__main__')",
  "        with log.open('a', encoding='utf-8') as stream:",
  "            stream.write('DONE: '+script+'\\n')",
  "    with log.open('a', encoding='utf-8') as stream:",
  "        stream.write('ANALYSES_COMPLETED\\n')",
  "    with zipfile.ZipFile(temporary_zip, 'w', zipfile.ZIP_DEFLATED) as archive:",
  "        for folder in ['results_water', 'results_climate']:",
  "            for file in sorted((root/folder).rglob('*')):",
  "                if file.is_file():",
  "                    archive.write(file, file.relative_to(root))",
  "        archive.write(log, log.name)",
  "    os.replace(temporary_zip, root/'RESULTADOS_PARA_DEVOLVER.zip')",
  "    with log.open('a', encoding='utf-8') as stream:",
  "        stream.write('ZIP_READY\\n')",
  "    print('\\nConcluido:', root/'RESULTADOS_PARA_DEVOLVER.zip', flush=True)",
  "except Exception:",
  "    with log.open('a', encoding='utf-8') as stream:",
  "        stream.write('FAILED\\n'+traceback.format_exc())",
  "    if temporary_zip.exists():",
  "        temporary_zip.unlink()",
  "    raise",
  "finally:",
  "    sys.argv = original_argv",
  sep = "\n"
)
reticulate::py_run_string(codigo_python, local = FALSE)

message(
  "Concluído. Envie este arquivo: ",
  file.path(revision_dir, "RESULTADOS_PARA_DEVOLVER.zip")
)
