from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from itertools import zip_longest
import re
from math import inf

try:
    import morphoTreeAdjust as mta  # package; exposes binding classes
except Exception:  # pragma: no cover
    mta = None  # for type hints; at runtime, install the [viz] extra


def _is_dynamic_tree(tree) -> bool:
    return all(hasattr(tree, attr) for attr in (
        "nodes",
        "root",
        "getChildren",
        "getAltitude",
        "getPixelsOfCC",
        "reconstructionImage",
    ))


def _as_image(image, *, num_rows=None, num_cols=None):
    array = np.asarray(image)
    if array.ndim == 1 and num_rows is not None and num_cols is not None:
        return array.reshape(num_rows, num_cols)
    return array


def _dynamic_tree_pixels(tree, node_id, *, use_subtree=False):
    if use_subtree and hasattr(tree, "getNodeSubtree"):
        pixels = []
        for subtree_node in tree.getNodeSubtree(node_id):
            pixels.extend(tree.getPixelsOfCC(subtree_node))
        return pixels
    return list(tree.getPixelsOfCC(node_id))


def _dynamic_tree_value(tree, node_id, attrs):
    values = []
    for attr in attrs:
        if attr == "id":
            values.append(node_id)
        elif attr in {"level", "altitude"}:
            values.append(tree.getAltitude(node_id))
        elif attr == "area":
            values.append(len(tree.getPixelsOfCC(node_id)))
        elif attr == "parent":
            values.append(tree.getNodeParent(node_id))
        elif attr == "children":
            values.append(tree.getChildren(node_id))
        elif attr == "num_children":
            values.append(len(tree.getChildren(node_id)))
        elif attr == "subtree_size" and hasattr(tree, "getNodeSubtree"):
            values.append(len(tree.getNodeSubtree(node_id)))
        else:
            values.append(f"<unsupported:{attr}>")
    return ": ".join(str(value) for value in values)



# --- Minimal ANSI-aware utilities ---
_ANSI_RE = re.compile(r"\x1b\[[0-?]*[ -/]*[@-~]")
_RESET = "\x1b[0m"


def _strip_ansi(text: str) -> str:
    return _ANSI_RE.sub("", text)


def _text_width(text: str) -> int:
    # Simplified width calculation (ignores ANSI sequences)
    return len(_strip_ansi(text))


def _ljust(text: str, amount: int, padding: str = ' ') -> str:
    pad = max(0, amount - _text_width(text))
    return text + (padding * pad)


class _NodeFormatter:
    """Formats a node as a text block with known width."""

    @classmethod
    def from_string(cls, content: str):
        lines = content.split('\n')
        width = max(_text_width(line) for line in lines)
        return cls(lines, width=width)

    def __init__(self, lines, *, width: int, middle_width: int = None):
        self.lines = lines
        self.width = width
        self.middle_width = middle_width

    def color_bg(self, color: str, add_space: bool) -> None:
        if add_space:
            self.lines = [f'{ color } { _ljust(line, self.width) } { _RESET }' for line in self.lines]
            self.width += 2
        else:
            self.lines = [f'{ color }{ _ljust(line, self.width) }{ _RESET }' for line in self.lines]

    def to_str(self) -> str:
        return '\n'.join(self.lines)

    def get_middle_width(self) -> int:
        if self.middle_width is None:
            return sum(divmod(self.width, 2)) - 1
        return self.middle_width


def _join_boxes(boxes):
    lines = [
        ' '.join(_ljust(line, boxes[i].width) for i, line in enumerate(lines))
        for lines in zip_longest(*(box.lines for box in boxes), fillvalue='')
    ]
    width = sum(box.width for box in boxes) + len(boxes) - 1
    return lines, width


def _add_pipes(boxes, lines) -> int:
    padding = ' ' * boxes[0].get_middle_width()
    pipes = '┌'
    for prev, box in zip(boxes, boxes[1:]):
        pipes += '─' * (prev.width - prev.get_middle_width() + box.get_middle_width()) + '┬'
    middle_of_pipes = sum(divmod(len(pipes), 2)) - 1
    pipes = (
        padding + pipes[:middle_of_pipes]
        + {"─": "┴", "┬": "┼", "┌": "├", "┐": "┤"}[pipes[middle_of_pipes]]
        + pipes[middle_of_pipes + 1:-1] + '┐'
    )
    lines.insert(0, pipes)
    return len(padding) + middle_of_pipes


def _join_horizontally(boxes):
    lines, width = _join_boxes(boxes)
    middle = _add_pipes(boxes, lines)
    return _NodeFormatter(lines, width=width, middle_width=middle)


def _add_parent(parent, children):
    parent_middle, children_middle = parent.get_middle_width(), children.get_middle_width()
    parent_width, children_width = parent.width, children.width
    if parent_middle == children_middle:
        lines = parent.lines + children.lines
        middle = parent_middle
    elif parent_middle < children_middle:
        padding = ' ' * (children_middle - parent_middle)
        lines = [padding + line for line in parent.lines] + children.lines
        parent_width += children_middle - parent_middle
        middle = children_middle
    else:
        padding = ' ' * (parent_middle - children_middle)
        lines = parent.lines + [padding + line for line in children.lines]
        children_width += parent_middle - children_middle
        middle = parent_middle
    return _NodeFormatter(lines, width=max(parent_width, children_width), middle_width=middle)


class _TreeFormatter:
    """Formats the tree into lines (vertical only)."""

    def __init__(self, get_children, get_val, color):
        self.get_children = get_children
        self.get_node_val = get_val
        self.color = color

    def format(self, node):
        return self._tree_vertical(node).to_str().rstrip()

    def _tree_vertical(self, node):
        children = self.get_children(node)
        cur = self._format_node(node)
        if children:
            children_fmt = [self._tree_vertical(child) for child in children]
            if len(children_fmt) == 1:
                only = children_fmt[0]
                only.lines.insert(0, ' ' * only.get_middle_width() + '|')
                children_node = only
            else:
                children_node = _join_horizontally(children_fmt)
            cur = _add_parent(cur, children_node)
        return cur

    def _format_node(self, node):
        contents = str(self.get_node_val(node))
        nf = _NodeFormatter.from_string(contents)
        if self.color:
            nf.color_bg(self.color, add_space=True)
        return nf


class PrintTree:
    """Prints the tree to the console in vertical format."""

    def __init__(self, get_children=None, get_val=None, *, color: str = ""):
        self.default_get_children = get_children or (lambda x: x.children)
        self.default_get_node_val = get_val or (lambda x: x.value)
        self.default_color = color

    def __call__(self, node, *, color: str = None):
        fmt = _TreeFormatter(
            get_children=self.default_get_children,
            get_val=self.default_get_node_val,
            color=self.default_color if color is None else color,
        )
        print(fmt.format(node))


def printTree(root, attrs=['id', 'level', 'area']):
    if _is_dynamic_tree(root):
        tree = root
        _printTree = PrintTree(
            lambda node_id: tree.getChildren(node_id),
            lambda node_id: _dynamic_tree_value(tree, node_id, attrs),
            color='\x1b[40m\x1b[37m'
        )
        return _printTree(tree.root)

    _printTree = PrintTree(
        lambda n: n.children,
        lambda n: ": ".join(str(getattr(n, a)) for a in attrs),
        color='\x1b[40m\x1b[37m'
    )
    return _printTree(root)


def computeCentroid(pixels, num_cols):
    n = len(pixels)
    
    # 1) continuous centroid
    sum_r = sum(p // num_cols for p in pixels)
    sum_c = sum(p %  num_cols for p in pixels)
    rc_mean = sum_r / n
    cc_mean = sum_c / n

    # 2) set pixel closest to the centroid
    best_idx = None
    best_r = best_c = None
    best_d2 = inf
    for p in pixels:
        r = p // num_cols
        c = p %  num_cols
        d2 = (r - rc_mean)**2 + (c - cc_mean)**2
        if d2 < best_d2 or (d2 == best_d2 and (best_idx is None or p < best_idx)):
            best_d2 = d2
            best_idx = p
            best_r, best_c = r, c

    return best_r, best_c

def showLevelSets(img_f):
    print("Upper and lower level sets where in the black (highlighted are cnps) are foreground pixels and white are background pixels")
    t_values = np.unique(img_f)
    len = np.size(t_values)
    max_value = np.max(img_f)
    fig, axes = plt.subplots(len, 2, figsize=(5, 15))

    # Fill the subplots with thresholded images
    for i, t in enumerate(t_values):
        # First column: img_f >= t
        thr = t_values[i]
        img_threshold_ge = np.where(img_f > thr, 1, 0)
        img_threshold_cnps = np.where(img_f == thr, 1, 0)

        axes[i, 0].imshow(img_threshold_ge, cmap='gray_r', vmax=1, vmin=0, interpolation='nearest')
        axes[i, 0].imshow(img_threshold_cnps, cmap='gray_r', vmax=1, vmin=0, alpha=0.4, interpolation='nearest')
        axes[i, 0].set_title(f'image ≥ {t}')
        axes[i, 0].axis('off')
        axes[i, 0].add_patch(plt.Rectangle((0, 0), img_threshold_ge.shape[1] - 1, img_threshold_ge.shape[0] - 1,
                                   edgecolor='red', linewidth=0.5, fill=False))
        # Second column: img_f <= t
        thr = t_values[len - i -1]
        img_threshold_le = np.where(img_f < thr, 1, 0)
        img_threshold_cnps = np.where(img_f == thr, 1, 0)
        axes[i, 1].imshow(img_threshold_le, cmap='gray_r',vmax=1, vmin=0, interpolation='nearest')
        axes[i, 1].imshow(img_threshold_cnps, cmap='gray_r', vmax=1, vmin=0, alpha=0.4, interpolation='nearest')
        axes[i, 1].set_title(f'image ≤ {thr}')
        axes[i, 1].axis('off')
        axes[i, 1].add_patch(plt.Rectangle((0, 0), img_threshold_le.shape[1] - 1, img_threshold_le.shape[0] - 1,
                                   edgecolor='blue', linewidth=0.5, fill=False))

    # Adjust layout to avoid overlap
    plt.tight_layout()
    plt.show()


def showTree(tree):
    if _is_dynamic_tree(tree):
        if(tree.isMaxtree):
            print("Image and its max-tree representation. The indices of max-tree nodes are shown as labels on the image.")
        else:
            print("Image and its min-tree representation. The indices of min-tree nodes are shown as labels on the image.")

        image = _as_image(tree.reconstructionImage(), num_rows=tree.numRows, num_cols=tree.numCols)
        ids_position = []
        for node_id in tree.nodes:
            pixels = tree.getPixelsOfCC(node_id)
            r, c = computeCentroid(pixels, tree.numCols)
            ids_position.append((node_id, r, c))

        plt.figure(figsize=(8, 6))
        ax = plt.gca()
        ax.add_patch(plt.Rectangle((0, 0), image.shape[1] - 1, image.shape[0] - 1, edgecolor='red', linewidth=0.5, fill=False))
        plt.imshow(image, cmap='gray')

        for node_id, y, x in ids_position:
            plt.text(x, y, f'•{node_id}', color='red', fontsize=9)

        plt.axis('off')
        plt.show()
        printTree(tree)
        return

    if(tree.isMaxtree):
        print("Image and its max-tree representation. The indices of max-tree nodes are shown as labels on the image.")
    else:
        print("Image and its min-tree representation. The indices of min-tree nodes are shown as labels on the image.")
    img_vector = tree.reconstructionImage()
    ids_position = []

    for node in tree.root.postOrderTraversal():
        for repFZ in node.repCNPs:
            fz = tree.getPixelsOfFlatzone(repFZ)
            r, c = computeCentroid(fz, tree.numCols)
            ids_position.append( (node.id, r, c) )

    # Plot the image
    plt.figure(figsize=(8, 6))
    image = img_vector.reshape(tree.numRows, tree.numCols)
    
    # Add a border around the image
    ax = plt.gca()  # Get the current axis
    ax.add_patch(plt.Rectangle((0, 0), image.shape[1] - 1, image.shape[0] - 1, edgecolor='red', linewidth=0.5, fill=False))
    plt.imshow(image, cmap='gray')

    # Add labels at the synthetic centroids
    for i, (id, y, x) in enumerate(ids_position):
        plt.text(x, y, f'•{id}', color='red', fontsize=9)


    plt.axis('off')
    plt.show()
    printTree(tree.root)

    
def showNode(tree, node, showCNPs=False):
    if _is_dynamic_tree(tree):
        pixels = _dynamic_tree_pixels(tree, node, use_subtree=False)
        image = _as_image(tree.reconstructionImage(), num_rows=tree.numRows, num_cols=tree.numCols)
        mask = np.zeros((tree.numRows, tree.numCols), dtype=np.uint8)
        for p in pixels:
            row, col = divmod(p, tree.numCols)
            mask[row, col] = 1

        fig, ax = plt.subplots(figsize=(6, 6))
        ax.imshow(image, cmap='gray_r')
        ax.imshow(mask, cmap='autumn', alpha=0.55)
        ax.add_patch(plt.Rectangle((0, 0), image.shape[1] - 1, image.shape[0] - 1, edgecolor='red', linewidth=0.5, fill=False))

        if showCNPs:
            for p in pixels:
                row, col = divmod(p, tree.numCols)
                ax.scatter(col, row, color='red', s=6)

        ax.set_title(
            f"node={node} level={tree.getAltitude(node)} "
            f"pixels={len(pixels)}"
        )
        ax.axis("off")
        plt.show()
        return

    image = tree.reconstructionNode(node).reshape(tree.numRows, tree.numCols)
    ax = plt.gca()  # Get the current axis
    ax.imshow(image, cmap='gray_r')
    ax.add_patch(plt.Rectangle((0, 0), image.shape[1] - 1, image.shape[0] - 1, edgecolor='red', linewidth=0.5, fill=False))
    if showCNPs:
        num_cols = tree.numCols
        for p in node.cnps:
            row, col = divmod(p, num_cols)  # Convert a linear index to 2D coordinates
            ax.scatter(col, row, color='red', s=5)


    ax.set_title(f"{node}")
    ax.axis("off")


def showSubtree(tree, node, showCNPs=False):
    if _is_dynamic_tree(tree):
        pixels = _dynamic_tree_pixels(tree, node, use_subtree=True)
        image = _as_image(tree.reconstructionImage(), num_rows=tree.numRows, num_cols=tree.numCols)
        mask = np.zeros((tree.numRows, tree.numCols), dtype=np.uint8)
        for p in pixels:
            row, col = divmod(p, tree.numCols)
            mask[row, col] = 1

        fig, ax = plt.subplots(figsize=(8, 6))
        ax.imshow(image, cmap='gray_r')
        ax.imshow(mask, cmap='autumn', alpha=0.55)
        ax.add_patch(plt.Rectangle((0, 0), image.shape[1] - 1, image.shape[0] - 1, edgecolor='red', linewidth=0.5, fill=False))

        if showCNPs:
            for p in pixels:
                row, col = divmod(p, tree.numCols)
                ax.scatter(col, row, color='red', s=6)

        subtree_size = len(tree.getNodeSubtree(node)) if hasattr(tree, "getNodeSubtree") else 1
        ax.set_title(
            f"subtree root={node} level={tree.getAltitude(node)} "
            f"subtree_size={subtree_size}"
        )
        ax.axis("off")
        plt.show()
        return

    return showNode(tree, node, showCNPs=showCNPs)


def treeRows(tree):
    if not _is_dynamic_tree(tree):
        raise TypeError("treeRows currently supports only the dynamic tree API.")

    rows = []
    for node_id in tree.nodes:
        row = {
            "node": node_id,
            "altitude": tree.getAltitude(node_id),
            "parent": tree.getNodeParent(node_id),
            "children": tree.getChildren(node_id),
            "num_pixels": len(tree.getPixelsOfCC(node_id)),
            "is_leaf": tree.isLeaf(node_id),
        }
        if hasattr(tree, "getNodeSubtree"):
            row["subtree_size"] = len(tree.getNodeSubtree(node_id))
        rows.append(row)
    return rows
