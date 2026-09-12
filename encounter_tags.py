"""Date-based encounter association; not an instrument operating-mode flag."""
SOURCE='https://science.nasa.gov/blogs/parker-solar-probe/2025/03/25/nasas-parker-solar-probe-completes-23rd-close-approach-to-sun/'
def tag_encounter(row):
    day=row['utc'][:10]
    # Reviewed coverage only: do not silently extrapolate to future inputs.
    if '2025-03-04' <= day <= '2025-03-31':
        period='before' if day<'2025-03-18' else 'after' if day>'2025-03-27' else 'within'
        label={'before':'E23 context — before encounter','within':'E23 — encounter dates','after':'E23 context — after encounter'}[period]
        row.update(psp_encounter=23,encounter_period=period,encounter_tag=label,encounter_perihelion_date='2025-03-22',encounter_source=SOURCE)
    else:
        row.update(psp_encounter=None,encounter_period='unassigned',encounter_tag='Encounter unassigned',encounter_perihelion_date=None,encounter_source=None)
    return row
