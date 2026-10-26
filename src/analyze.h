#ifndef ANALYZE_H
#define ANALYZE_H
/* C6 driver: conjunction events CSV (from screen/fastscreen) -> graph, clusters, Top-K, mitigation.
 * movable_csv: OMM CSV whose NORAD ids are the manoeuvrable objects (e.g. data/raw/active.csv), or NULL = everything movable.
 * Writes <prefix>_components.csv, <prefix>_topk.csv, <prefix>_cover.csv and prints a summary plus a machine-readable
 * "C6SUMMARY,..." line. Returns 0, or non-zero on error. */
int analyze_events(const char *events_csv, const char *movable_csv, int K, const char *prefix, long exact_node_limit);
#endif
