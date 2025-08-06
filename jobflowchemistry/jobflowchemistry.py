"""
Primary functions for the JobFlowChemistry package.
"""


def canvas(with_attribution=True):
    """
    Return a motivational quote about code as a placeholder example.

    Parameters
    ----------
    with_attribution : bool, optional
        Whether to include the attribution (default is True).

    Returns
    -------
    str
        The quote, optionally with attribution.
    """
    quote = "The code is but a canvas to our imagination."
    if with_attribution:
        quote += "\n\t- Adapted from Henry David Thoreau"
    return quote


if __name__ == "__main__":
    # Do something if this file is invoked on its own
    print(canvas())
