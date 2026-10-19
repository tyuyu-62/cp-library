# Library: https://github.com/tyuyu-62/cp-library
_HASH_MOD = (1 << 61) - 1


def rooted_tree_isomorphism_classification_safe(
    parent: list[int],
) -> tuple[int, list[int]]:
    """Classify rooted subtrees exactly.

    The parent array may describe a rooted forest. Roots use -1 and every
    other vertex must have a smaller parent index. The number of vertices
    must be at most 500000.

    Time: O(N + sum(degree(v) log degree(v))).
    Space: O(N).
    """
    n = len(parent)
    if n == 0:
        return 0, []
    child_count = bytearray(n)
    first = [0] * n
    second = [0] * n
    extra = [None] * n
    classification = [0] * n
    table = {}
    intern = table.setdefault
    pair_tag = 1 << 38
    unary = [0] * n
    kinds = 1
    for v in range(n - 1, -1, -1):
        count = child_count[v]
        if count == 0:
            kind = 0
        elif count == 1:
            child_kind = first[v]
            kind = unary[child_kind]
            if kind == 0:
                kind = kinds
                unary[child_kind] = kind
                kinds += 1
        elif count == 2:
            key = pair_tag | (first[v] << 19) | second[v]
            kind = intern(key, kinds)
            if kind == kinds:
                kinds += 1
        else:
            signature = extra[v]
            signature.sort()
            key = tuple(signature)
            extra[v] = None
            kind = intern(key, kinds)
            if kind == kinds:
                kinds += 1
        classification[v] = kind
        p = parent[v]
        if p >= 0:
            count = child_count[p]
            if count == 0:
                first[p] = kind
                child_count[p] = 1
            elif count == 1:
                previous = first[p]
                if previous <= kind:
                    second[p] = kind
                else:
                    first[p] = kind
                    second[p] = previous
                child_count[p] = 2
            elif count == 2:
                extra[p] = [first[p], second[p], kind]
                child_count[p] = 3
            else:
                extra[p].append(kind)
    return kinds, classification


def rooted_tree_isomorphism_classification_fast(
    parent: list[int],
) -> tuple[int, list[int]]:
    """Classify rooted subtrees with a randomized polynomial fingerprint.

    The parent array may describe a rooted forest. Roots use -1 and every
    other vertex must have a smaller parent index. The number of vertices
    must be at most 500000.

    Time: expected O(N).
    Space: O(N).

    Paths and trees with one branching vertex whose children are leaves are
    classified directly. Trees of maximum degree at most two use the safe
    classifier. Other trees use the product over children of (child hash + a
    random weight indexed by child subtree size), modulo 2**61 - 1. A fresh
    256-bit secret seeds the pseudorandom generator on each call. Assuming
    its outputs act as independent uniform field elements, the union bound
    for any collision in a fixed forest of at most 500000 vertices is below
    5.7e-7. The bound uses that equal-sized rooted subtrees are disjoint.
    Use the safe function when a deterministic guarantee is required.
    """
    n = len(parent)
    if n == 0:
        return 0, []
    for v in range(1, n):
        if parent[v] != v - 1:
            break
    else:
        return n, list(range(n))
    if n < 1024:
        return rooted_tree_isomorphism_classification_safe(parent)
    if parent[0] < 0:
        child_count = bytearray(n)
        first_child = [-1] * n
        branch = -1
        multiple_branches = False
        nonleaf_branch_child = False
        degree_over_two = False
        root_count = 0
        for v in range(n - 1, -1, -1):
            p = parent[v]
            if p < 0:
                root_count += 1
                continue
            if branch >= 0:
                if p == branch and child_count[v]:
                    nonleaf_branch_child = True
                elif parent[p] == branch:
                    nonleaf_branch_child = True
            count = child_count[p]
            if count == 0:
                first_child[p] = v
            elif count == 1:
                if branch < 0:
                    branch = p
                    if child_count[first_child[p]] or child_count[v]:
                        nonleaf_branch_child = True
                elif branch != p:
                    multiple_branches = True
            elif count == 2:
                degree_over_two = True
            if count < 3:
                child_count[p] = count + 1
            if degree_over_two and (
                multiple_branches or nonleaf_branch_child
            ):
                break
        if (
            branch >= 0
            and not multiple_branches
            and not nonleaf_branch_child
            and root_count == 1
        ):
            classification = [0] * n
            v = branch
            kind = 1
            while v >= 0:
                classification[v] = kind
                kind += 1
                v = parent[v]
            return kind, classification
        if not degree_over_two:
            return rooted_tree_isomorphism_classification_safe(parent)
    import random
    import secrets

    mod = _HASH_MOD
    sizes = [1] * n
    hashes = [1] * n
    weights = [None] * (n + 1)
    rng = random.Random(secrets.randbits(256))
    multiply = _mul_mod
    for v in range(n - 1, -1, -1):
        p = parent[v]
        if p >= 0:
            size = sizes[v]
            sizes[p] += size
            weight = weights[size]
            if weight is None:
                weight = rng.getrandbits(61)
                while weight == mod:
                    weight = rng.getrandbits(61)
                weights[size] = weight
            value = hashes[v] + weight
            if value >= mod:
                value -= mod
            hashes[p] = multiply(hashes[p], value)
    table = {}
    intern = table.get
    classification = [0] * n
    kinds = 0
    for v in range(n):
        key = (sizes[v], hashes[v])
        kind = intern(key)
        if kind is None:
            kind = kinds
            table[key] = kind
            kinds += 1
        classification[v] = kind
    return kinds, classification


def _mul_mod(a: int, b: int) -> int:
    mod = _HASH_MOD
    au, ad = a >> 31, a & 0x7FFFFFFF
    bu, bd = b >> 31, b & 0x7FFFFFFF
    m = ad * bu + au * bd
    mu, md = m >> 30, m & 0x3FFFFFFF
    x = (au * bu << 1) + mu + (md << 31) + ad * bd
    xu, xd = x >> 61, x & mod
    result = xu + xd
    return result if result < mod else result - mod
