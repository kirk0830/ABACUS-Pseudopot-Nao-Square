from reportlab.pdfgen import canvas
from datetime import datetime
from pytz import timezone

def title(c: canvas.Canvas) -> int:
    """
    Create title for the pdf file, return the lowest y coordinate of all the elements in the title
    """
    # add title: "ABACUS Pseudopotential Square", centered
    c.setFont("Helvetica", 20)
    c.drawCentredString(300, 700, "ABACUS Pseudopotential Square")
    # add image "logo_basic.png" in current folder at (50, 600) and scale to (100, 100)
    c.drawInlineImage("logo_basic.png", 100, 700, 50, 50)
    # add text "Generated automatically based on ReportLab", italic, gray, at bottom right corner of title
    c.setFont("Helvetica-Oblique", 10)
    c.setFillColorRGB(0.5, 0.5, 0.5)
    c.drawRightString(550, 680, "Generated automatically based on Python.ReportLab")
    # add time and time zone information below the string above

    c.setFont("Helvetica", 10)
    c.setFillColorRGB(0, 0, 0)
    fmt = "%Y-%m-%d %H:%M:%S %Z%z"
    now_time = datetime.now(timezone("Asia/Shanghai"))
    c.drawRightString(550, 660, now_time.strftime(fmt))

    return 600

def paragraph(c: canvas.Canvas, y: int, title: str, content: str) -> int:
    """
    Add a paragraph of text to the pdf file, starting from y, return the lowest y coordinate of all the elements in the paragraph
    """
    # add title, bold, left aligned, starting from y
    c.setFont("Helvetica-Bold", 12)
    c.drawString(50, y, title)
    # add content, left aligned, starting from y-20
    c.setFont("Helvetica", 10)
    c.drawString(50, y-20, content)
    # return y-40
    return y-40

def main():
    c = canvas.Canvas("hello.pdf")
    y = title(c)
    content = """
To create a figure with two x-axes, you can use the `twinx()` method of the `Axes` object. Here's an example:

```python
import matplotlib.pyplot as plt

# Create some dummy data
x1 = [1, 2, 3, 4, 5]
y1 = [10, 8, 6, 4, 2]
x2 = [2, 4, 6, 8, 10]
y2 = [5, 7, 9, 11, 13]

# Create the figure and the first axis
fig, ax1 = plt.subplots()

# Plot the first data set on the first axis
ax1.plot(x1, y1, 'b-')
ax1.set_xlabel('X-axis 1 label')
ax1.set_ylabel('Y-axis label 1', color='b')

# Create the second axis and plot the second data set on it
ax2 = ax1.twinx()
ax2.plot(x2, y2, 'r-')
ax2.set_xlabel('X-axis 2 label')
ax2.set_ylabel('Y-axis label 2', color='r')

# Show the plot
plt.show()
```

In this example, we create a figure with one axis (`ax1`) and plot the first data set on it. We then create a second axis (`ax2`) using the `twinx()` method and plot the second data set on it. We set the labels and colors for each axis separately using the `set_xlabel()`, `set_ylabel()`, and `color` arguments. Finally, we show the plot using `plt.show()`.
"""
    y = paragraph(c, y, "Introduction", content)

    c.save()

if __name__ == "__main__":
    main()