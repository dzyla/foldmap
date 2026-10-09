"""Place sheet blocks and helices on the page (page units: 1 unit = one strand pitch)."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .features import helix_bundles
from .frame import Frame
from .model import SSE, Nucleic, Sheet, duplex_axis, duplex_groups
from .style import Style

PITCH = 1.1  # distance between neighbouring strand centres
SHAFT = 0.62  # arrow shaft width
SCALE = 0.23  # page units per Å (PITCH / 4.8 Å strand spacing)
HEAD = 1.0  # strand arrowhead width
HELIX_W = 1.0  # helix ribbon diameter, ~4.4 Å: one arrowhead wide, as in a cartoon (Style.helix_scale)
STRAND_RISE = 3.3  # Å per residue along a strand
HELIX_RISE = 1.5  # Å per residue along a helix
TEN_RISE = 2.0  # Å per residue along a 3-10 helix
ETA_W = 0.7  # 3-10 box width as a fraction of the helix width
MARGIN = 0.8  # clear space kept around each element so loops can pass
_ROW_GAP = 2.5  # between neighbours in a stack row
_ROW_SPACING = 3.2  # between stack rows: room for a chain-end label below one row and above the next
LABEL_GAP = 0.3  # a helix label's centre sits this far beyond the side of the element box
END_STUB = 0.7  # chain-end stub length
END_LABEL = 1.15  # distance from a chain-end port to its N/C letter
_GLYPH_W, _GLYPH_H = 0.32, 0.42  # rough label glyph size in page units, for keep-out boxes
_MAX_ROW = 30.0
_MIN_STRAND = 1.0  # shortest drawn strand (room for the arrow head)
_W_SEQ = 1.0  # pull between the two ends of a loop
_LOOP_SLACK = 1.5  # loop ends closer than this are left alone
_W_SS = 2.0  # pull between two bridged cysteines (twice a loop: disulfides should read as short bars)
_SS_SLACK = 0.6
_W_CONTACT = 0.3  # pull between elements that touch in 3D (at _CONTACT_CAP CA pairs)
_CONTACT_CAP = 20
_K_HOME = 0.02  # pull back toward the projected position, keeps the fold recognisable
_STEP = 0.15
_SWEEPS = 400
_SEP_SWEEPS = 5000  # push-apart passes before the guaranteed slide
_LAMBDAS = (0.3, 3.0, 30.0, 300.0)  # overlap penalty, stage by stage
_MAX_MOVE = 0.5  # largest step an item takes in one sweep


@dataclass(eq=False)
class Placed:
    sse: SSE
    cx: float
    cy: float
    length: float
    angle: float  # direction N -> C, radians (0 = right, pi/2 = up)
    width: float
    label: str = ""
    ghost: bool = False  # repeat of a barrel's first strand, drawn to show the sheet closes
    ends: str = ""  # "N", "C" or "NC" when the chain starts/ends here (its stub and letter need room)

    @property
    def direction(self) -> np.ndarray:
        return np.array([np.cos(self.angle), np.sin(self.angle)])

    @property
    def n_port(self) -> tuple[float, float]:
        p = np.array([self.cx, self.cy]) - self.direction * self.length / 2
        return float(p[0]), float(p[1])

    @property
    def c_port(self) -> tuple[float, float]:
        p = np.array([self.cx, self.cy]) + self.direction * self.length / 2
        return float(p[0]), float(p[1])

    @property
    def exit_c(self) -> np.ndarray:
        return self.direction

    @property
    def exit_n(self) -> np.ndarray:
        return -self.direction

    @property
    def rect(self) -> tuple[float, float, float, float]:
        hx = abs(np.cos(self.angle)) * self.length / 2 + abs(np.sin(self.angle)) * self.width / 2
        hy = abs(np.sin(self.angle)) * self.length / 2 + abs(np.cos(self.angle)) * self.width / 2
        return (self.cx - hx, self.cy - hy, self.cx + hx, self.cy + hy)


DNA_RISE = 3.4  # Å per base pair (B-DNA)
DNA_W = 20.0 * SCALE  # B-DNA diameter
DNA_PITCH = 10.5 * DNA_RISE * SCALE  # one turn, page units
DNA_MINOR = 0.37  # the second backbone trails the first by this fraction of a turn (minor groove)


@dataclass(eq=False)
class PlacedDNA:
    """A duplex (or lone strand) drawn as a straight double helix. axial[(strand, nt)] is the nucleotide's
    position along the drawn axis from the centre, in page units (positive = along `direction`)."""

    id: str
    strands: list[int]
    cx: float
    cy: float
    length: float
    angle: float
    width: float
    axial: dict[tuple[int, int], float]
    ghost: bool = False

    @property
    def direction(self) -> np.ndarray:
        return np.array([np.cos(self.angle), np.sin(self.angle)])

    rect = Placed.rect


@dataclass(eq=False)
class Layout:
    placed: dict[str, Placed]
    bounds: tuple[float, float, float, float]
    ghosts: list[Placed] = field(default_factory=list)
    dna: list[PlacedDNA] = field(default_factory=list)
    nucleic: Nucleic | None = None
    sheets: list[list[str]] = field(default_factory=list)  # strand ids of each sheet block, in layout order
    domains: list[tuple[str, list[str]]] = field(default_factory=list)  # named domains and their element ids
    nt_rows: bool = False  # base letters printed as sequence rows above each duplex
    focus: set[str] | None = None  # chains drawn in colour when highlighting; the rest in grey
    links: object = None  # model.Links: disulfides and glycans, residue indices in Backbone order
    res_chain: list[str] = field(default_factory=list)  # chain of each Backbone residue
    res_name: list[str] = field(default_factory=list)  # residue type of each Backbone residue
    res_b: list[float] = field(default_factory=list)  # CA B-factor of each Backbone residue
    membrane: dict | None = None  # {"y": (bottom, top) of the band, "inside": "below"|"above", "source": ...}
    res_cons: list[float] = field(default_factory=list)  # alignment conservation of each residue (nan: none)

    @property
    def width(self) -> float:
        return self.bounds[2] - self.bounds[0]

    @property
    def height(self) -> float:
        return self.bounds[3] - self.bounds[1]


@dataclass(eq=False)
class _Item:
    """A sheet block or a single helix, moved as a unit. members hold offsets from (cx, cy)."""

    members: list[Placed]
    offsets: list[tuple[float, float]]
    cx: float
    cy: float
    home: tuple[float, float]
    order: int  # first residue index, for stack mode
    chain: str  # chain of the lowest-numbered element, for stack mode

    @property
    def box(self) -> tuple[float, float, float, float]:
        """Footprint relative to (cx, cy): every member with its label room."""
        fp = [_footprint(m) for m in self.members]
        return (
            min(o[0] + f[0] for o, f in zip(self.offsets, fp)),
            min(o[1] + f[1] for o, f in zip(self.offsets, fp)),
            max(o[0] + f[2] for o, f in zip(self.offsets, fp)),
            max(o[1] + f[3] for o, f in zip(self.offsets, fp)),
        )

    @property
    def half(self) -> tuple[float, float]:
        x0, y0, x1, y1 = self.box
        return (x1 - x0) / 2, (y1 - y0) / 2

    @property
    def mid(self) -> tuple[float, float]:
        """Centre of the footprint relative to (cx, cy)."""
        x0, y0, x1, y1 = self.box
        return (x0 + x1) / 2, (y0 + y1) / 2

    def apply(self) -> None:
        for m, (ox, oy) in zip(self.members, self.offsets):
            m.cx, m.cy = self.cx + ox, self.cy + oy


def _footprint(m) -> tuple[float, float, float, float]:
    """A member's box relative to its own centre; a helix's includes the room its label needs."""
    x0, y0, x1, y1 = m.rect
    box = (x0 - m.cx, y0 - m.cy, x1 - m.cx, y1 - m.cy)
    if isinstance(m, Placed) and m.sse.kind in ("H", "G"):
        lx, ly = helix_label_pos(Placed(m.sse, 0.0, 0.0, m.length, m.angle, m.width))
        hx, hy = _GLYPH_W * 1.5 + 0.05, _GLYPH_H / 2  # up to three glyphs: α10
        box = (min(box[0], lx - hx), min(box[1], ly - hy), max(box[2], lx + hx), max(box[3], ly + hy))
    if isinstance(m, Placed) and m.ends:
        d = m.direction
        for end in m.ends:
            out = d if end == "C" else -d
            t = out * (m.length / 2 + END_LABEL)  # the N/C letter beyond the stub
            hx, hy = _GLYPH_W / 2 + 0.05, _GLYPH_H / 2
            box = (min(box[0], t[0] - hx), min(box[1], t[1] - hy), max(box[2], t[0] + hx), max(box[3], t[1] + hy))
    return box


def _strand_length(sse: SSE) -> float:
    return max(_MIN_STRAND, (len(sse) - 1) * STRAND_RISE * SCALE)


def _anchor(strands: list[SSE]) -> tuple[int, str]:
    first = min(strands, key=lambda s: s.start)
    return first.start, first.chain


def _sheet_item(sheet: Sheet, frame: Frame, pitch: float = PITCH) -> _Item:
    strands = sheet.strands
    cents = frame.project(np.array([s.centroid for s in strands]), [s.chain for s in strands])
    mean_dir = sum(d * s.axis for d, s in zip(sheet.directions, strands))
    at = np.mean([x.centroid for x in strands], axis=0)
    seen = frame.direction(at, mean_dir)
    s = 1 if seen[1] >= 0 else -1  # +1: strands[0] points up
    # Drawing the strands vertical turns the projected sheet in the page plane; the strand order must turn
    # with it (left to right = across the drawn-up direction), or the sheet shows its back face (mirror image).
    up2d = s * seen
    right = np.array([up2d[1], -up2d[0]]) if np.linalg.norm(up2d) > 1e-6 else np.array([1.0, 0.0])
    order = list(range(len(strands)))
    if len(strands) > 1 and (cents[-1] - cents[0]) @ right < 0:
        order.reverse()
    members, offsets = [], []
    slots = len(order) + (1 if sheet.closed else 0)
    for slot, k in enumerate(order):
        up = sheet.directions[k] * s > 0
        sse = strands[k]
        members.append(Placed(sse, 0, 0, _strand_length(sse), np.pi / 2 if up else -np.pi / 2, pitch))
        offsets.append(((slot - (slots - 1) / 2) * pitch, 0.0))
    if sheet.closed:  # the barrel wraps: its first strand pairs again with the last
        f = members[0]
        members.append(Placed(f.sse, 0, 0, f.length, f.angle, f.width, ghost=True))
        offsets.append(((slots - 1) / 2 * pitch, 0.0))
    home = tuple(cents.mean(axis=0) * SCALE)
    return _Item(members, offsets, home[0], home[1], home, *_anchor(strands))


def _helix_item(
    sse: SSE, frame: Frame, angle_mode: str = "snap", width: float = HELIX_W, rise: float = HELIX_RISE
) -> _Item:
    """Always true length. angle_mode: 'tilted' as seen in the view; 'snap' to the nearest of up, down, left,
    right (the topology-diagram convention); 'upright' up or down like a strand (stack mode)."""
    au, av = frame.direction(sse.centroid, sse.axis)
    length = max(width, (len(sse) - 1) * rise * SCALE)  # true length even when tilted: equal pitch everywhere
    seen = float(np.arctan2(av, au))
    if angle_mode == "upright":
        angle = np.pi / 2 if av >= 0 else -np.pi / 2
    elif angle_mode == "snap":
        angle = float(np.round(seen / (np.pi / 2)) * (np.pi / 2))
    else:
        angle = seen
    pos = frame.project(sse.centroid, [sse.chain])[0] * SCALE
    placed = Placed(sse, 0, 0, length, angle, width)
    return _Item(
        [placed], [(0.0, 0.0)], float(pos[0]), float(pos[1]), (float(pos[0]), float(pos[1])), sse.start, sse.chain
    )


def _angle(seen: float, av: float, mode: str) -> float:
    if mode == "upright":
        return np.pi / 2 if av >= 0 else -np.pi / 2
    if mode == "snap":
        return float(np.round(seen / (np.pi / 2)) * (np.pi / 2))
    return seen


def _bundle_item(group: list[SSE], frame: Frame, angle_mode: str, width: float) -> _Item:
    """Helices of several chains packed along their length, drawn as one unit: one column per helix (kinked
    pieces of one helix share a column, end to end), columns side by side in their order across the bundle at
    their true spacing (at least a label's width apart), staggered along the shared axis as in 3D;
    antiparallel helices point the other way."""
    from .features import KINK_GAP

    ref = group[0].axis
    mean = sum(np.sign(float(s.axis @ ref)) * s.axis * len(s) for s in group)
    m = mean / np.linalg.norm(mean)
    centre = np.mean([s.centroid for s in group], axis=0)
    a2 = frame.direction(centre, m)
    angle = _angle(float(np.arctan2(a2[1], a2[0])), float(a2[1]), angle_mode)
    d2, side2 = np.array([np.cos(angle), np.sin(angle)]), np.array([-np.sin(angle), np.cos(angle)])
    columns: list[list[SSE]] = []
    for s in sorted(group, key=lambda s: (s.chain, s.start)):
        last = columns[-1][-1] if columns else None
        if last is not None and last.chain == s.chain and s.start - last.end - 1 <= KINK_GAP:
            columns[-1].append(s)  # the next piece of a kinked helix
        else:
            columns.append([s])
    lateral = [float(np.mean([frame.project(s.centroid, [s.chain])[0] @ side2 for s in col])) for col in columns]
    columns = [c for _, c in sorted(zip(lateral, columns), key=lambda x: x[0])]
    pts = [np.mean([s.centroid for s in col], axis=0) for col in columns]
    gaps = [np.linalg.norm((a - b) - ((a - b) @ m) * m) for a, b in zip(pts, pts[1:])]
    spacing = max(width + 1.2, float(np.mean(gaps)) * SCALE if gaps else 0.0)
    members, offsets = [], []
    for k, col in enumerate(columns):
        last_end = -np.inf
        for s in sorted(col, key=lambda s: float((s.centroid - centre) @ m)):
            length = max(width, (len(s) - 1) * HELIX_RISE * SCALE)
            along = max(float((s.centroid - centre) @ m) * SCALE, last_end + 0.3 + length / 2)  # kinked pieces in a row
            last_end = along + length / 2
            members.append(Placed(s, 0, 0, length, angle if float(s.axis @ m) >= 0 else angle + np.pi, width))
            lat = (k - (len(columns) - 1) / 2) * spacing
            offsets.append(tuple(float(v) for v in d2 * along + side2 * lat))
    pos = frame.project(centre)[0] * SCALE
    first = min(group, key=lambda s: s.start)
    return _Item(
        members, offsets, float(pos[0]), float(pos[1]), (float(pos[0]), float(pos[1])), first.start, first.chain
    )


def _separate(items: list[_Item]) -> None:
    """Push overlapping items apart (axis of least penetration); spring back to home early on."""
    halves = [it.half for it in items]
    mids = [it.mid for it in items]
    for sweep in range(_SEP_SWEEPS):
        moved = False
        for i in range(len(items)):
            for j in range(i + 1, len(items)):
                a, b = items[i], items[j]
                dx, dy = b.cx + mids[j][0] - a.cx - mids[i][0], b.cy + mids[j][1] - a.cy - mids[i][1]
                px = halves[i][0] + halves[j][0] + MARGIN - abs(dx)
                py = halves[i][1] + halves[j][1] + MARGIN - abs(dy)
                if px > 1e-9 and py > 1e-9:
                    moved = True
                    if px < py:
                        s = (1 if dx > 0 else -1) if dx != 0 else 1
                        a.cx -= s * px / 2
                        b.cx += s * px / 2
                    else:
                        s = (1 if dy > 0 else -1) if dy != 0 else 1
                        a.cy -= s * py / 2
                        b.cy += s * py / 2
        if sweep < 100:
            for it in items:
                it.cx += 0.03 * (it.home[0] - it.cx)
                it.cy += 0.03 * (it.home[1] - it.cy)
        if not moved and sweep >= 100:
            return
    # guarantee: slide remaining collisions to the right, one item at a time
    placed: list[int] = []
    for i in range(len(items)):
        while any(_collide(items[i], halves[i], mids[i], items[j], halves[j], mids[j]) for j in placed):
            items[i].cx += 0.5
        placed.append(i)


def _owner(items: list[_Item]) -> dict[str, tuple[int, Placed, tuple[float, float]]]:
    return {
        (m.id if isinstance(m, PlacedDNA) else m.sse.id): (k, m, o)
        for k, it in enumerate(items)
        for m, o in zip(it.members, it.offsets)
        if not m.ghost
    }


def _contact_weights(items: list[_Item], contacts: dict[frozenset[str], int]) -> np.ndarray:
    """Symmetric pull strength between items whose elements touch in 3D."""
    owner = _owner(items)
    w = np.zeros((len(items), len(items)))
    for pair, count in sorted((tuple(sorted(k)), v) for k, v in contacts.items()):
        a, b = (owner.get(x) for x in pair)
        if a and b and a[0] != b[0]:
            w[a[0], b[0]] += _W_CONTACT * min(count, _CONTACT_CAP) / _CONTACT_CAP
            w[b[0], a[0]] = w[a[0], b[0]]
    return w


def _loop_springs(items: list[_Item], sses: list[SSE]) -> list[tuple[int, np.ndarray, int, np.ndarray]]:
    """(item of a, C-port offset of a, item of b, N-port offset of b) for each loop between items."""
    owner = _owner(items)
    springs = []
    for a, b in zip(sses, sses[1:]):
        if a.chain != b.chain or a.id not in owner or b.id not in owner:
            continue
        (i, ma, oa), (j, mb, ob) = owner[a.id], owner[b.id]
        if i != j:
            springs.append((i, np.add(oa, ma.direction * ma.length / 2), j, np.add(ob, -mb.direction * mb.length / 2)))
    return springs


def _mass(it: _Item) -> float:
    """How firmly an item holds its projected place: one per element; a duplex like a sheet of
    one strand per two base pairs (the protein moves around the DNA, not the DNA around the protein)."""
    return sum(len(m.axial) / 4 if isinstance(m, PlacedDNA) else 1 for m in it.members)


def _energy(flat, home, anchor_w, want, W, loops, lam, mids=None):
    """Layout energy and its gradient. Terms: stay near the projected position (fold); loop ends within
    _LOOP_SLACK of each other; touching elements' (and DNA binders') boxes adjacent;
    lam * squared box overlap (raised stage by stage until overlaps are gone)."""
    P = flat.reshape(-1, 2)
    n = len(P)
    grad = 2 * anchor_w * (P - home)
    E = float((anchor_w[:, 0] * ((P - home) ** 2).sum(1)).sum())
    Q = P if mids is None else P + mids  # footprint centres
    d = Q[None, :, :] - Q[:, None, :]  # d[i, j] = Q_j - Q_i
    sg = np.sign(d)
    px, py = want[..., 0] - np.abs(d[..., 0]), want[..., 1] - np.abs(d[..., 1])  # overlap depth per axis
    upper = np.triu(np.ones((n, n), bool), 1)

    def pair_grad(gx, gy):  # dE/d(overlap depth) per pair (upper triangle) -> per item
        grad[:, 0] += (gx * sg[..., 0]).sum(1) - (gx * sg[..., 0]).sum(0)
        grad[:, 1] += (gy * sg[..., 1]).sum(1) - (gy * sg[..., 1]).sum(0)

    on_x = px < py
    ov = upper & (px > 0) & (py > 0)
    depth = np.where(on_x, px, py)
    E += lam * float((depth[ov] ** 2).sum())
    pair_grad(np.where(ov & on_x, 2 * lam * px, 0.0), np.where(ov & ~on_x, 2 * lam * py, 0.0))
    gap = np.maximum(-px, -py)  # how far apart two boxes are beyond touching
    ct = upper & (W > 0) & (gap > 0)
    E += float((W[ct] * gap[ct] ** 2).sum())
    cg = np.where(ct, 2 * W * gap, 0.0)
    gx_side = -px >= -py
    pair_grad(-np.where(gx_side, cg, 0.0), -np.where(~gx_side, cg, 0.0))
    I, OI, J, OJ, w, slack = loops
    if len(I):
        g = (P[J] + OJ) - (P[I] + OI)
        dist = np.maximum(np.hypot(g[:, 0], g[:, 1]), 1e-12)
        excess = np.clip(dist - slack, 0.0, None)
        E += float((w * excess**2).sum())
        dg = (2 * w * excess / dist)[:, None] * g
        np.add.at(grad, J, dg)
        np.add.at(grad, I, -dg)
    return E, grad.ravel()


def _attract(
    items: list[_Item],
    sses: list[SSE],
    contacts: dict[frozenset[str], int],
    ties: dict[int, tuple[int, int]] | None = None,
    spread: float = 0.0,
    lift: float = 0.0,
    bridges: list[tuple[str, float, str, float]] = (),
) -> None:
    """Minimise the layout energy (_energy) from the projected positions, raising the overlap penalty in
    stages. The result becomes each item's new home; _separate then guarantees no overlap remains."""
    from scipy.optimize import minimize

    n = len(items)
    W = _contact_weights(items, contacts)
    springs = _loop_springs(items, sses)
    weight = [_W_SEQ] * len(springs)
    slack = [_LOOP_SLACK] * len(springs)
    owner = _owner(items)
    bridges_used = []
    for a, ta, b, tb in bridges:  # disulfides: the two cysteines themselves pulled together
        if a in owner and b in owner and owner[a][0] != owner[b][0]:
            bridges_used.append((a, b))
            (i, ma, oa), (j, mb, ob) = owner[a], owner[b]
            springs.append(
                (
                    i,
                    np.add(oa, ma.direction * (ta - 0.5) * ma.length),
                    j,
                    np.add(ob, mb.direction * (tb - 0.5) * mb.length),
                )
            )
            weight.append(_W_SS)
            slack.append(_SS_SLACK)
    loops = (
        np.array([s[0] for s in springs], int),
        np.array([s[1] for s in springs], float).reshape(-1, 2),
        np.array([s[2] for s in springs], int),
        np.array([s[3] for s in springs], float).reshape(-1, 2),
        np.array(weight, float),
        np.array(slack, float),
    )
    home = np.array([it.home for it in items], float)
    h = np.array([it.half for it in items], float)
    want = h[:, None, :] + h[None, :, :] + MARGIN
    anchor_w = _K_HOME * np.array([[_mass(it)] for it in items], float)
    mids = np.array([it.mid for it in items], float)
    ties = ties or {}
    free = [i for i in range(n) if i not in ties]
    slot = {i: k for k, i in enumerate(free)}
    base = np.array([slot[ties[i][0]] if i in ties else slot[i] for i in range(n)])
    copy_k = np.array([ties[i][1] if i in ties else 0 for i in range(n)], float)

    def expand(z):
        P = z[:-1].reshape(-1, 2)[base].copy()
        P[:, 0] += copy_k * z[-1]  # protomer k sits k spreads to the right of protomer 0
        P[:, 1] += copy_k * lift  # and, along a filament, k rises higher
        return P

    def tied(z):  # one position per orbit of symmetry copies, plus the spread between protomers
        E, g = _energy(expand(z).ravel(), home, anchor_w, want, W, loops, lam, mids)
        g = g.reshape(-1, 2)
        gz = np.zeros((len(free), 2))
        np.add.at(gz, base, g)
        return E, np.append(gz.ravel(), float((g[:, 0] * copy_k).sum()))

    def solve(z, stages=_LAMBDAS):
        nonlocal lam
        for lam in stages:  # noqa: B007 - read by tied() through nonlocal
            z = minimize(tied, z, jac=True, method="L-BFGS-B", options={"maxiter": 300}).x
        return z

    lam = _LAMBDAS[0]
    z = solve(np.append(np.array([[items[i].cx, items[i].cy] for i in free], float).ravel(), spread))
    lam = _LAMBDAS[-1]
    best = tied(z)[0]
    mass = np.array([_mass(it) for it in items])
    for k in range(len(springs) - len(bridges_used), len(springs)):  # disulfides still stretched: a second start
        i, oi, j, oj = springs[k]
        P = expand(z)
        if np.hypot(*((P[j] + oj) - (P[i] + oi))) <= 2 * _SS_SLACK + 1.0:
            continue
        mover, anchor, om, oa = (i, j, oi, oj) if mass[i] <= mass[j] else (j, i, oj, oi)
        if mover in ties:
            continue
        trial = z.copy()
        target = P[anchor] + oa - om + np.array([0.0, -1.0])  # the lighter element's cysteine beside its partner
        trial[2 * slot[mover] : 2 * slot[mover] + 2] = target
        trial = solve(trial, _LAMBDAS[1:])
        lam = _LAMBDAS[-1]
        energy = tied(trial)[0]
        if energy < best:
            z, best = trial, energy
    x = expand(z).ravel()
    for it, p in zip(items, x.reshape(n, 2)):
        it.cx, it.cy = float(p[0]), float(p[1])
        it.home = (it.cx, it.cy)


def _collide(a: _Item, ha, ma, b: _Item, hb, mb) -> bool:
    """Do the footprints (half extents h, centre offsets m) of a and b, plus MARGIN, overlap?"""
    dx, dy = b.cx + mb[0] - a.cx - ma[0], b.cy + mb[1] - a.cy - ma[1]
    return ha[0] + hb[0] + MARGIN - abs(dx) > 1e-9 and ha[1] + hb[1] + MARGIN - abs(dy) > 1e-9


def _stack(items: list[_Item]) -> None:
    """Rows in sequence order, one row break per chain and when a row is full."""
    items.sort(key=lambda it: (it.chain, it.order))
    rows: list[list[_Item]] = []
    width = 0.0
    for it in items:
        hx = it.half[0]
        new_chain = rows and rows[-1][-1].chain != it.chain
        if not rows or new_chain or width + 2 * hx > _MAX_ROW:
            rows.append([])
            width = 0.0
        rows[-1].append(it)
        width += 2 * hx + _ROW_GAP
    y_top = 0.0
    for row in rows:
        height = max(it.box[3] - it.box[1] for it in row)
        x = 0.0
        for it in row:  # footprints (with label room) side by side, centred in the row's height
            it.cx = x - it.box[0]
            it.cy = y_top - height / 2 - it.mid[1]
            x += it.box[2] - it.box[0] + _ROW_GAP
        dna_row = all(isinstance(m, PlacedDNA) for it in row for m in it.members)
        y_top -= height + (MARGIN if dna_row else _ROW_SPACING)  # binders hang right under the DNA


def _centre_dna_over_binders(items: list[_Item], contacts: dict[frozenset[str], int]) -> None:
    owner = _owner(items)
    for it in items:
        dna = [m for m in it.members if isinstance(m, PlacedDNA)]
        if not dna:
            continue
        xs = [
            items[owner[x][0]].cx + owner[x][2][0]
            for key in contacts
            if dna[0].id in key
            for x in key
            if x != dna[0].id and x in owner
        ]
        if xs:
            it.cx = float(np.mean(xs))


def _assign_labels(placed: list[Placed]) -> None:
    letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    counts: dict[tuple[str, str], int] = {}
    for p in sorted(placed, key=lambda p: (p.sse.chain, p.sse.start)):
        n = counts.get((p.sse.chain, p.sse.kind), 0)
        counts[(p.sse.chain, p.sse.kind)] = n + 1
        if p.sse.kind == "E":
            p.label = letters[n] if n < 26 else f"S{n + 1}"
        elif p.sse.kind == "G":
            p.label = f"η{n + 1}"
        else:
            p.label = f"α{n + 1}"


def _duplexes(na: Nucleic) -> list[list[int]]:
    return duplex_groups(na)


def _dna_item(k: int, strands: list[int], na: Nucleic, frame: Frame, flat: bool, width: float = DNA_W) -> _Item:
    centre, axis = duplex_axis(na, strands)  # the first strand runs 5' -> 3' along the drawn direction
    along = {
        (s, i): float((na.strands[s].c1[i] - centre) @ axis) * SCALE  # C1' sits level across a base pair
        for s in strands
        for i in range(len(na.strands[s]))
    }
    for a, i, b, j in na.pairs:  # a base pair is drawn as one straight rung
        if (a, i) in along and (b, j) in along:
            along[(a, i)] = along[(b, j)] = (along[(a, i)] + along[(b, j)]) / 2
    span = max(along.values()) - min(along.values())
    mid = (max(along.values()) + min(along.values())) / 2
    along = {key: v - mid for key, v in along.items()}
    length = span + DNA_RISE * SCALE  # half a step beyond the end base pairs on each side
    a2 = frame.direction(centre, axis)
    angle = 0.0 if flat or np.linalg.norm(a2) < 1e-6 else float(np.arctan2(a2[1], a2[0]))
    pos = frame.project(centre)[0] * SCALE
    d = PlacedDNA(f"dna:{k}", list(strands), 0.0, 0.0, length, angle, width, along)
    return _Item([d], [(0.0, 0.0)], float(pos[0]), float(pos[1]), (float(pos[0]), float(pos[1])), -1, "")


def _dna_contacts(sses: list[SSE], na: Nucleic, groups: list[list[int]]) -> dict[frozenset[str], int]:
    """Residues of each element touching each duplex, as extra adjacency pull between them."""
    out: dict[frozenset[str], int] = {}
    for k, strands in enumerate(groups):
        for s in sses:
            n = sum(1 for r in range(s.start, s.end + 1) if any(x in strands for x, _ in na.contacts.get(r, ())))
            if n:
                out[frozenset((s.id, f"dna:{k}"))] = 4 * n  # one residue in the groove counts like several CA pairs
    return out


def _symmetry_ties(items: list[_Item], frame) -> dict[int, tuple[int, int]]:
    """Copy item -> (its protomer-0 mate, protomer index k): matched by chain class and shared residue numbers.
    The copy also takes the mate's shape and orientation, so every protomer is drawn the same."""
    sector, classes = getattr(frame, "sector", None), getattr(frame, "classes", None)
    if not sector or not classes:
        return {}

    def residues(it):  # chains outside every protomer (bound peptides, unpaired chains) take no part
        return {
            (classes[m.sse.chain], r)
            for m in it.members
            if isinstance(m, Placed) and not m.ghost and m.sse.chain in classes
            for r in range(m.sse.first.seq, m.sse.last.seq + 1)
        }

    def protomer(it):
        ks = {sector.get(m.sse.chain) for m in it.members if isinstance(m, Placed)}
        return ks.pop() if len(ks) == 1 else None

    sets = [residues(it) for it in items]
    zero = [i for i, it in enumerate(items) if protomer(it) == 0 and sets[i]]
    ties = {}
    for i, it in enumerate(items):
        k = protomer(it)
        if not k or not sets[i]:
            continue
        best = max(zero, key=lambda j: len(sets[i] & sets[j]) / len(sets[i] | sets[j]), default=None)
        if best is None or len(sets[i] & sets[best]) < 0.5 * len(sets[i] | sets[best]):
            continue
        mate = items[best]

        def span(m):
            if m.sse.chain not in classes:
                return set()
            return {(classes[m.sse.chain], r) for r in range(m.sse.first.seq, m.sse.last.seq + 1)}

        pairs = [
            max(
                range(len(mate.members)),
                key=lambda q: (len(span(m) & span(mate.members[q])), mate.members[q].ghost == m.ghost),
            )
            for m in it.members
        ]
        if len(mate.members) == len(it.members) and sorted(pairs) == list(range(len(pairs))):
            for m, q in zip(it.members, pairs):  # one-to-one: take the mate's shape member by member
                m.angle, m.length = mate.members[q].angle, mate.members[q].length
            it.offsets = [mate.offsets[q] for q in pairs]
        else:  # DSSP split or merged an element in one copy: keep the shape, match the orientation
            votes = [
                np.cos(m.angle - mate.members[q].angle)
                for m, q in zip(it.members, pairs)
                if span(m) & span(mate.members[q])
            ]
            if votes and np.mean(votes) < 0:  # drawn upside down: turn the block 180 degrees in the page
                for m in it.members:
                    m.angle += np.pi
                it.offsets = [(-ox, -oy) for ox, oy in it.offsets]
        ties[i] = (best, k)
    return ties


def resolve(layout: Layout, sses: list[SSE], ref: str) -> list[str]:
    """Element ids for a reference: a label (G, α2; every chain), chain:label (A:G), #k (the k-th element
    from the N terminus, 1-based), or res:a-b / res:chain:a-b (elements overlapping those residue numbers)."""
    ref = ref.strip()
    placed = [s for s in sses if s.id in layout.placed]
    if ref.startswith("#"):
        k = int(ref[1:]) if ref[1:].isdigit() else 0
        if not 1 <= k <= len(placed):
            raise ValueError(f"no element {ref!r}: elements are numbered #1 to #{len(placed)}")
        return [placed[k - 1].id]
    if ref.startswith("res:"):
        body = ref[4:]
        chain = None
        if body.count(":") == 1:
            chain, body = body.split(":")
        a, _, b = body.partition("-")
        try:
            lo, hi = int(a), int(b or a)
        except ValueError:
            raise ValueError(f"residue reference {ref!r} must look like res:150-159 or res:A:150-159") from None
        out = [s.id for s in placed if (chain is None or s.chain == chain) and s.first.seq <= hi and s.last.seq >= lo]
    elif ":" in ref:
        chain, label = ref.split(":", 1)
        out = [s.id for s in placed if s.chain == chain and layout.placed[s.id].label == label]
    else:
        out = [s.id for s in placed if layout.placed[s.id].label == ref]
    if not out:
        raise ValueError(f"no element {ref!r}; use a label such as A or α1, chain:label, #k or res:a-b")
    return out


def _one(layout: Layout, sses: list[SSE], ref: str) -> str:
    ids = resolve(layout, sses, ref)
    if len(ids) != 1:
        raise ValueError(f"{ref!r} matches {len(ids)} elements ({', '.join(ids)}); name one, e.g. chain:label")
    return ids[0]


def _adjust(
    items: list[_Item],
    lay: Layout,
    sses: list[SSE],
    rename: dict[str, str],
    swap: list[tuple[str, str]],
    move: list[tuple[str, tuple[float, float]]],
) -> None:
    """The user's own hand on the layout: new labels, exchanged places, nudges; then overlaps are cleared."""
    where = {
        m.sse.id: (it, k) for it in items for k, m in enumerate(it.members) if isinstance(m, Placed) and not m.ghost
    }
    for ref, name in rename.items():
        for sid in resolve(lay, sses, ref):
            it, k = where[sid]
            it.members[k].label = name
            for g in it.members:
                if g.ghost and g.sse.id == sid:
                    g.label = name
    moved = False
    for a_ref, b_ref in swap:
        a, b = _one(lay, sses, a_ref), _one(lay, sses, b_ref)
        (ia, ka), (ib, kb) = where[a], where[b]
        if ia is ib:  # two strands of one sheet: exchange their slots
            ia.offsets[ka], ia.offsets[kb] = ia.offsets[kb], ia.offsets[ka]
        else:
            ia.cx, ib.cx, ia.cy, ib.cy = ib.cx, ia.cx, ib.cy, ia.cy
        moved = True
    for ref, (dx, dy) in move:
        for sid in resolve(lay, sses, ref):
            it, _ = where[sid]
            it.cx += dx
            it.cy += dy
            moved = True
    if moved:
        for it in items:
            it.home = (it.cx, it.cy)
        _separate(items)
        for it in items:
            it.apply()


def build_layout(
    sses: list[SSE],
    sheets: list[Sheet],
    frame: Frame,
    mode: str = "projected",
    *,
    attract: bool = True,
    contacts: dict[frozenset[str], int] | None = None,
    nucleic: Nucleic | None = None,
    style: Style | None = None,
    bridges: list[tuple[str, float, str, float]] | None = None,
    domains: list[tuple[str, list[str]]] | None = None,
    rename: dict[str, str] | None = None,
    swap: list[tuple[str, str]] | None = None,
    move: list[tuple[str, tuple[float, float]]] | None = None,
) -> Layout:
    """`attract` (projected mode) pulls loop-joined and touching elements together before overlaps
    are removed; `contacts` comes from features.sse_contacts. `nucleic` adds DNA/RNA duplexes, placed
    like any element and pulling the elements that touch them."""
    if mode not in ("projected", "stack"):
        raise ValueError(f"unknown layout mode {mode!r}")
    style = style or Style()
    angle_mode = "upright" if mode == "stack" else style.helix_angle
    items = [_sheet_item(sh, frame, PITCH * style.strand_scale) for sh in sheets]
    bundles = helix_bundles(sses)
    bundled = {k for g in bundles for k in g}
    by_id = {s.id: s for s in sses}
    items += [_bundle_item([by_id[k] for k in g], frame, angle_mode, HELIX_W * style.helix_scale) for g in bundles]
    items += [
        _helix_item(s, frame, angle_mode, HELIX_W * style.helix_scale)
        for s in sses
        if s.kind == "H" and s.id not in bundled
    ]
    items += [
        _helix_item(s, frame, angle_mode, HELIX_W * ETA_W * style.helix_scale, TEN_RISE) for s in sses if s.kind == "G"
    ]
    if not items:
        return Layout({}, (0.0, 0.0, 0.0, 0.0))
    contacts = dict(contacts or {})
    if nucleic is not None and nucleic.strands:
        groups = _duplexes(nucleic)
        items += [
            _dna_item(k, g, nucleic, frame, mode == "stack", DNA_W * style.dna_scale) for k, g in enumerate(groups)
        ]
    if nucleic is not None and nucleic.strands:
        for key, n in _dna_contacts(sses, nucleic, groups).items():
            contacts[key] = contacts.get(key, 0) + n
    items.sort(key=lambda it: it.order)
    members = {m.sse.id: m for it in items for m in it.members if isinstance(m, Placed) and not m.ghost}
    for chain in dict.fromkeys(s.chain for s in sses):
        mine = [s for s in sses if s.chain == chain and s.id in members]
        if mine:
            members[mine[0].id].ends += "N"
            members[mine[-1].id].ends += "C"
    if mode == "stack":
        _stack(items)
        _centre_dna_over_binders(items, contacts)
    else:
        ties = _symmetry_ties(items, frame)
        spread, lift = 0.0, 0.0
        if ties:  # protomers side by side at least one protomer's width apart; copies start on their mates
            mates = {j for j, _ in ties.values()}
            x0 = min(items[j].cx + items[j].box[0] for j in mates)
            x1 = max(items[j].cx + items[j].box[2] for j in mates)
            spread = max(frame.period * SCALE, x1 - x0 + 2 * MARGIN)
            lift = getattr(frame, "rise", 0.0) * SCALE
            for i, (j, k) in ties.items():
                items[i].cx, items[i].cy = items[j].cx + k * spread, items[j].cy + k * lift
                items[i].home = (items[j].home[0] + k * spread, items[j].home[1] + k * lift)
        if attract and len(items) > 1:
            _attract(items, sses, contacts, ties, spread, getattr(frame, "rise", 0.0) * SCALE, bridges or [])
        _separate(items)
        if domains and _part_domains(items, domains):
            _separate(items)
    for it in items:
        it.apply()
    dna = [m for it in items for m in it.members if isinstance(m, PlacedDNA)]
    placed = [m for it in items for m in it.members if isinstance(m, Placed) and not m.ghost]
    ghosts = [m for it in items for m in it.members if isinstance(m, Placed) and m.ghost]
    _assign_labels(placed)
    real = {p.sse.id: p for p in placed}
    if rename or swap or move:
        _adjust(items, Layout(real, (0.0, 0.0, 0.0, 0.0)), sses, rename or {}, swap or [], move or [])
    for g in ghosts:
        g.label = real[g.sse.id].label
    rects = np.array([p.rect for p in placed + ghosts + dna])
    sheet_blocks = [
        [m.sse.id for m in it.members if not m.ghost]
        for it in items
        if it.members and all(isinstance(m, Placed) and m.sse.kind == "E" for m in it.members)
    ]
    lay = Layout(
        real,
        (0.0, 0.0, 0.0, 0.0),
        ghosts,
        dna,
        nucleic if dna else None,
        sheet_blocks,
        bool(dna) and style.nucleotide_labels,
    )
    lay.domains = [(name, [k for k in ids if k in real]) for name, ids in domains or []]
    rects = np.vstack([rects, *[[domain_panel(lay, name)] for name, _ in lay.domains]]) if lay.domains else rects
    rects = np.vstack([rects, [b for _, b in label_boxes(lay, sses)] or np.empty((0, 4))])
    lay.bounds = (
        float(rects[:, 0].min()),
        float(rects[:, 1].min()),
        float(rects[:, 2].max()),
        float(rects[:, 3].max()),
    )
    return lay


DOMAIN_PAD = 0.7  # panel margin around a domain's elements
_DOMAIN_GAP = 2 * DOMAIN_PAD + 0.6  # clear space between two domains' elements: both panels plus a gutter


def _part_domains(items: list[_Item], domains: list[tuple[str, list[str]]], sweeps: int = 60) -> bool:
    """Move whole domains apart (each as a rigid group) until their panels no longer overlap. Each item goes
    with the domain most of its elements belong to. True if anything moved."""
    owner = {k: d for d, (_, ids) in enumerate(domains) for k in ids}
    groups: dict[int, list[_Item]] = {}
    for it in items:
        votes = [owner[m.sse.id] for m in it.members if isinstance(m, Placed) and m.sse.id in owner]
        if votes:
            groups.setdefault(max(set(votes), key=votes.count), []).append(it)
    keys = sorted(groups)
    if len(keys) < 2:
        return False

    def box(d):
        b = np.array([(it.cx + it.box[0], it.cy + it.box[1], it.cx + it.box[2], it.cy + it.box[3]) for it in groups[d]])
        return b[:, 0].min(), b[:, 1].min(), b[:, 2].max(), b[:, 3].max()

    moved = False
    for _ in range(sweeps):
        clash = False
        for i, a in enumerate(keys):
            for b in keys[i + 1 :]:
                A, B = box(a), box(b)
                ox = min(A[2], B[2]) - max(A[0], B[0]) + _DOMAIN_GAP
                oy = min(A[3], B[3]) - max(A[1], B[1]) + _DOMAIN_GAP
                if ox <= 0 or oy <= 0:
                    continue
                clash = moved = True
                ca, cb = (
                    np.array([(A[0] + A[2]) / 2, (A[1] + A[3]) / 2]),
                    np.array([(B[0] + B[2]) / 2, (B[1] + B[3]) / 2]),
                )
                if ox <= oy:  # the cheaper direction, each side going half way
                    step = np.array([ox / 2 if cb[0] >= ca[0] else -ox / 2, 0.0])
                else:
                    step = np.array([0.0, oy / 2 if cb[1] >= ca[1] else -oy / 2])
                for it, sign in [(it, -1) for it in groups[a]] + [(it, 1) for it in groups[b]]:
                    it.cx, it.cy = it.cx + sign * step[0], it.cy + sign * step[1]
                    it.home = (it.home[0] + sign * step[0], it.home[1] + sign * step[1])  # rest here from now on
        if not clash:
            break
    return moved


def domain_panel(layout: Layout, name: str) -> tuple[float, float, float, float]:
    """A domain's panel (its elements plus margin) with room above for its name."""
    ids = dict(layout.domains)[name]
    r = np.array([layout.placed[k].rect for k in ids] + [g.rect for g in layout.ghosts if g.sse.id in ids])
    x0, y0, x1, y1 = r[:, 0].min(), r[:, 1].min(), r[:, 2].max(), r[:, 3].max()
    return (float(x0 - DOMAIN_PAD), float(y0 - DOMAIN_PAD), float(x1 + DOMAIN_PAD), float(y1 + DOMAIN_PAD + 0.6))


def provisional(sses: list[SSE]) -> Layout:
    """A layout holding only element labels, for resolving references before anything is placed."""
    placed = [Placed(s, 0.0, 0.0, 1.0, 0.0, 1.0) for s in sses]
    _assign_labels(placed)
    return Layout({p.sse.id: p for p in placed}, (0.0, 0.0, 0.0, 0.0))


def helix_label_pos(p: Placed) -> tuple[float, float]:
    """Beside the bar, on its upper side, clear of the body."""
    nx, ny = -np.sin(p.angle), np.cos(p.angle)
    side = (p.width / 2 + LABEL_GAP) * (1 if ny >= 0 else -1)
    return float(p.cx + nx * side), float(p.cy + ny * side)


def termini(layout: Layout, sses: list[SSE]) -> list[tuple[str, str, tuple[float, float], np.ndarray]]:
    """(N or C, chain, port, outward direction) for each chain's first and last placed element."""
    out = []
    for chain in dict.fromkeys(s.chain for s in sses):
        mine = [s for s in sses if s.chain == chain and s.id in layout.placed]
        if mine:
            first, last = layout.placed[mine[0].id], layout.placed[mine[-1].id]
            out += [("N", chain, first.n_port, first.exit_n), ("C", chain, last.c_port, last.exit_c)]
    return out


def _box(cx: float, cy: float, hx: float, hy: float) -> tuple[float, float, float, float]:
    return (cx - hx, cy - hy, cx + hx, cy + hy)


def dna_y(slot: int, x: np.ndarray | float, width: float = DNA_W) -> np.ndarray | float:
    """Sideways offset of backbone `slot` (0 or 1) at axial position x on a drawn duplex (page units)."""
    return width / 2 * np.sin(2 * np.pi * (np.asarray(x) / DNA_PITCH - slot * DNA_MINOR))


def dna_to_page(d: PlacedDNA, x, y) -> np.ndarray:
    """Duplex-local (axial, sideways) to page coordinates."""
    c, s = np.cos(d.angle), np.sin(d.angle)
    x, y = np.asarray(x, float), np.asarray(y, float)
    return np.stack([d.cx + x * c - y * s, d.cy + x * s + y * c], axis=-1)


NT_ROW = 0.55  # page units between the two sequence rows (and from the duplex to the first)


def dna_letter_pos(d: PlacedDNA, row: int, x: float) -> tuple[float, float]:
    """Where a base letter goes: sequence rows just above the duplex, the first strand on top, its partner
    under it, one column per base pair."""
    p = dna_to_page(d, x, d.width / 2 + NT_ROW * (2 - row))
    return float(p[0]), float(p[1])


def dna_end_labels(d: PlacedDNA, na: Nucleic) -> list[tuple[str, str, tuple[float, float]]]:
    """(strand id, 5′ or 3′, page position) just beyond each strand's terminal nucleotides."""
    out = []
    reach = DNA_RISE * SCALE / 2 + 0.5
    for slot, k in enumerate(d.strands):
        n = len(na.strands[k])
        sign = 1 if slot % 2 == 0 else -1  # slot 0 runs 5' -> 3' along the drawn direction, slot 1 back
        for end, i in (("5′", 0), ("3′", n - 1)):
            x = d.axial[(k, i)] + (reach if (end == "3′") == (sign > 0) else -reach)
            pos = dna_to_page(d, x, dna_y(slot % 2, x, d.width))
            out.append((na.strands[k].id, end, (float(pos[0]), float(pos[1]))))
    return out


def label_boxes(layout: Layout, sses: list[SSE]) -> list[tuple[str, tuple[float, float, float, float]]]:
    """Keep-out boxes for text and chain-end stubs: loops route around them, the page bounds include them."""
    out = []
    for p in layout.placed.values():
        if p.sse.kind in ("H", "G"):
            out.append(
                (f"label:{p.sse.id}", _box(*helix_label_pos(p), _GLYPH_W * len(p.label) / 2 + 0.05, _GLYPH_H / 2))
            )
    for end, chain, port, ex in termini(layout, sses):
        a, b = np.add(port, ex * 0.15), np.add(port, ex * END_STUB)
        out.append(
            (
                f"{end}-stub:{chain}",
                (min(a[0], b[0]) - 0.05, min(a[1], b[1]) - 0.05, max(a[0], b[0]) + 0.05, max(a[1], b[1]) + 0.05),
            )
        )
        t = np.add(port, ex * END_LABEL)
        out.append((f"{end}:{chain}", _box(float(t[0]), float(t[1]), _GLYPH_W / 2 + 0.05, _GLYPH_H / 2)))
    for d in layout.dna:
        if layout.nt_rows:
            for row, k in enumerate(d.strands[:2]):
                xs = [d.axial[(k, i)] for i in range(len(layout.nucleic.strands[k]))]
                (x0, y0), (x1, y1) = dna_letter_pos(d, row, min(xs)), dna_letter_pos(d, row, max(xs))
                out.append(
                    (
                        f"dna-seq:{k}",
                        (
                            min(x0, x1) - _GLYPH_W / 2,
                            min(y0, y1) - _GLYPH_H / 2,
                            max(x0, x1) + _GLYPH_W / 2,
                            max(y0, y1) + _GLYPH_H / 2,
                        ),
                    )
                )
        for sid, end, (x, y) in dna_end_labels(d, layout.nucleic):
            out.append((f"dna-end:{sid}:{end}", _box(x, y, _GLYPH_W + 0.05, _GLYPH_H / 2)))
    for name, ids in layout.domains:  # a domain's name, top left inside its panel
        if ids:
            x0, _, _, y1 = domain_panel(layout, name)
            out.append(
                (
                    f"domain-label:{name}",
                    (x0 + 0.25, y1 - 0.35 - _GLYPH_H * 1.1, x0 + 0.45 + _GLYPH_W * 1.1 * len(name), y1 - 0.2),
                )
            )
    return out
