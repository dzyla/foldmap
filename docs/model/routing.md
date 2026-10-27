# Loop routing

## Scope

- `route.route_loops`: the loop between each pair of sequence-adjacent elements, routed around everything placed.
  When the layout holds nucleic acid it also routes the single strands (`route_na.route_nucleic`, see
  [nucleic-acids.md](nucleic-acids.md)).
- `route._Router`: orthogonal A* on a grid of `route.RES` (0.25 page units), each drawing checked against the
  elements.
- `route._clear_helix_names`: after routing, a helix name a line runs through is moved.
- `route._walled_fallback`: in a membrane view, a fallback curve that would cross its loop's wall.
- `fallback.note`: each fallback taken is counted per figure (`fallback.FALLBACKS`).

## Assumes

- Elements are placed first. Routing never moves an element.
- A loop leaves its element's C-terminal end and enters the next one's N-terminal end.
- Each search has a node budget: `route._MAX_NODES` (200,000). A search that runs out with its ends joined by free
  cells runs once more with `route._MAX_NODES_WIDE` (1,600,000). All second searches of one figure share
  `route._MAX_NODES_PAGE` (8,000,000); the fine windows of a coarse page have their own pool,
  `route._MAX_NODES_WINDOWS` (8,000,000).
- A loop with no route is drawn as a smooth curve and counted (`loop_curve`). A route found only by opening the
  keep-outs around text is drawn over that text and counted (`loop_over_text`).

## Guarantees

- Two strands are joined by an orthogonal path around obstacles:
  `tests/test_route.py::test_two_strands_are_joined_by_an_orthogonal_path_around_the_obstacles`.
- No route runs through an element in random layouts:
  `tests/test_route.py::test_random_layouts_never_route_through_an_element`.
- A route through an element that cannot be redrawn is marked a fallback:
  `tests/test_route.py::test_a_route_through_an_element_that_cannot_be_redrawn_is_marked_fallback`.
- A walled-in port falls back and finishes: `tests/test_route.py::test_walled_in_port_falls_back_and_finishes`.
- Routing is deterministic: `tests/test_route.py::test_deterministic`.
- The second searches share one budget per page, and a search out of budget spends exactly its budget:
  `tests/test_route.py::test_the_second_searches_share_one_budget_per_page`,
  `tests/test_route.py::test_a_search_out_of_budget_spends_exactly_its_budget`.
- Every route constant measured was stable at x0.75 and x1.25 on the third sample (`route._SLAB_RUN`,
  `route._MAX_NODES_WIDE`, `route._MAX_NODES_PAGE`, `route._MAX_NODES_WINDOWS`, `route._NA_TEXT_TOLL`,
  `route._REROUTE_TURNS_INK`; `docs/validation/sensitivity/report.md`). Stable means no key moved on that sample at
  those two values.
- On the three samples' default views, no second search was refused for the page's budget
  (`docs/validation/readability-floor/report.md`, census).

## Fails

- Route fallbacks (`A3`) on the third sample at 7696947: 8XJV (four nucleic-acid strands drawn as plain curves), 6QIK
  (strand y:8-9), 7F5S (domains view, strand L5:694-695), 9FZL (domains view, loops of 4E and 4K and their copies)
  (`docs/validation/random100-third/review.md`).
- 9FZL's domains view: 98 chains, the page held at its cap, 30 second searches refused for the budget. With the
  windows given their own pool (`route._MAX_NODES_WINDOWS`), 2 fallback curves remain and the window pool is spent.
- Loops over a label (`A4`): 9NMN (helix α153's name), 9FZL's domains view (four loops)
  (`docs/validation/random100-third/review.md`).
- 7NZ2's domains view (second sample) spends the page budget. Without `route._MAX_NODES_PAGE` it newly fails `A4` and
  `A6`. The constant met the prune condition on the measured sets and is held for this reason
  (`docs/validation/prune/report.md`).
- `route._MAX_NODES_WINDOWS` is set for 7NZ2 and 9FZL; no value of it is measured against drawing from the page's
  budget again (`docs/validation/prune/report.md`, "needs code removal").

## Seen by

- `A3`: route fallbacks (`loop_curve`, the nucleic-acid strand fallbacks, `loop_wall_detour`, `loop_wall_unkept`).
- `A4`: a loop over a label (`loop_over_text`).
- `A18`: in a membrane view, a loop against the band and its wall.
- `R1`: the page's route budget spent (`route_second_search_refused`).
- Not seen by any check: `route_unsnapped` (a route redrawn unsnapped because the snapped one crossed an element). It
  is counted; the inventory's check column says "none".
- Not seen by any check: `helix_name_moved` (a helix name moved off a line). It is counted, and no criterion gates it.
- Not seen by any check: loops crossing loops. `figcheck.loop_crossings` and `figcheck.honest_crossings` count them
  and the campaign reports them; no criterion judges them.
- Time: a route search that is slow but finishes is judged only by `A2` (the default view's wall time, campaign
  only), which moved between two runs of one code on one machine (`docs/validation/sensitivity/report.md`).

## Evidence

- [inventory/inventory.md](../validation/inventory/inventory.md): the named fallbacks and the check that reports each.
- [sensitivity/report.md](../validation/sensitivity/report.md): the route constants, all stable.
- [prune/report.md](../validation/prune/report.md): `route:_MAX_NODES_PAGE` held; `route:_MAX_NODES_WINDOWS` not tried.
- [readability-floor/report.md](../validation/readability-floor/report.md): the route budget census.
- [random100-third/review.md](../validation/random100-third/review.md): `A3` and `A4` failures.
