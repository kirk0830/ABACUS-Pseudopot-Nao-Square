def concatenate(figures: list, direction: str, remove_after_quit: bool = False):
    """concatenate figures in the given direction, horizontally or vertically.
    Figures are given as .png, .svg or .jpg, or .bmp files."""
    d = "v" if direction.lower() in ["v", "vertical"] else "h"
    import os
    from PIL import Image

    figures = [f for f in figures if os.path.isfile(f)]
    # sort figures by names
    figures.sort()
    images = [Image.open(f) for f in figures]
    widths, heights = zip(*(i.size for i in images))

    if d == "v":
        max_width = max(widths)
        total_height = sum(heights)
        new_im = Image.new('RGB', (max_width, total_height))
        y_offset = 0
        for im in images:
            new_im.paste(im, (0, y_offset))
            y_offset += im.size[1]

    elif d == "h":
        total_width = sum(widths)
        max_height = max(heights)
        new_im = Image.new('RGB', (total_width, max_height))
        x_offset = 0
        for im in images:
            new_im.paste(im, (x_offset, 0))
            x_offset += im.size[0]

    feos = f"eos_concatenated_{d}.png"
    new_im.save(feos)
    
    if remove_after_quit:
        for f in figures:
            os.remove(f)

    return feos