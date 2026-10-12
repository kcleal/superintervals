from array import array

from superintervals import IntervalMap


def test_search_values() -> None:
    imap = IntervalMap()
    imap.add(10, 20, "A")
    imap.add(15, 25, "B")
    imap.add(30, 40, "C")
    imap.build()
    assert sorted(imap.search_values(8, 20)) == ["A", "B"]
    assert imap.count(8, 20) == 2
    assert imap.search_values(26, 29) == []


def test_from_arrays_and_batch_queries() -> None:
    imap = IntervalMap.from_arrays(array("i", [10, 15, 30]), array("i", [20, 25, 40]), ["A", "B", "C"])
    assert imap.count_batch(array("i", [5, 18, 35]), array("i", [12, 22, 45])) == [1, 2, 1]


def test_has_overlaps_with_nested_intervals() -> None:
    imap = IntervalMap()
    imap.add(0, 99, "long")
    imap.add(10, 19, "nested")
    imap.build()
    assert imap.has_overlaps(30, 39)
    assert imap.has_overlaps(15, 15)
    assert not imap.has_overlaps(100, 200)
