"""Combine the six parameter figures as vector artwork on landscape A4 pages.
Run with the bundled Python (requires pypdf and reportlab). By default read and
write alongside this script; --directory and --output override those paths.
"""
from pathlib import Path
from io import BytesIO
import argparse
from pypdf import PdfReader, PdfWriter, Transformation
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import A4, landscape

FIGURES = [
 ('Strategy coexistence - beta planes', 'Figure_strategy_coexistence_fixed_vs_responsive.pdf'),
 ('Strategy coexistence - phi planes', 'Figure_strategy_coexistence_phi_space_by_beta_fixed_vs_responsive.pdf'),
 ('Invested energy - beta planes', 'Figure_energy_allocation_fixed_vs_responsive.pdf'),
 ('Invested energy - phi planes', 'Figure_energy_allocation_phispace_fixed_vs_responsive.pdf'),
 ('Equal-talent 50:50 crossings - beta planes', 'Figure_equal_talent_5050_crossings.pdf'),
 ('Equal-talent 50:50 crossings - phi planes', 'Figure_equal_talent_5050_crossings_phispace.pdf'),
]

def combine(directory, output):
    writer = PdfWriter()
    width, height = landscape(A4)
    for number, (title, filename) in enumerate(FIGURES, 1):
        reader = PdfReader(directory / filename)
        if len(reader.pages) != 1:
            raise ValueError(f'Expected one figure page in {filename}')
        figure = reader.pages[0]
        figure.transfer_rotation_to_content()
        box = figure.mediabox
        fw, fh = float(box.width), float(box.height)
        left_margin, bottom_margin, top_margin = 18, 17, 35
        available_width = width - 2 * left_margin
        available_height = height - bottom_margin - top_margin
        scale = min(available_width / fw, available_height / fh)
        x = (width - scale * fw) / 2
        y = bottom_margin + (available_height - scale * fh) / 2
        frame = BytesIO()
        c = canvas.Canvas(frame, pagesize=(width, height))
        c.setFont('Helvetica', 10)
        c.drawString(left_margin, height - 18, f'{number}. {title}')
        c.setFont('Helvetica', 8)
        c.setFillGray(.4)
        c.drawRightString(width - left_margin, 8, f'{number} / {len(FIGURES)}')
        c.showPage()
        c.save()
        frame.seek(0)
        page = PdfReader(frame).pages[0]
        transform = (Transformation().translate(-float(box.left), -float(box.bottom))
                         .scale(scale).translate(x, y))
        page.merge_transformed_page(figure, transform, over=True)
        writer.add_page(page)
        writer.add_outline_item(title, number - 1)
    writer.add_metadata({'/Title': 'Metatournament parameter figures - print collection',
                         '/Subject': 'Six figures, one per landscape A4 page',
                         '/Creator': 'combine_parameter_figures_for_print.py'})
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open('wb') as stream:
        writer.write(stream)
    check = PdfReader(output)
    assert len(check.pages) == len(FIGURES)
    assert all(float(p.mediabox.width) > float(p.mediabox.height) for p in check.pages)
    print(f'Saved {output} ({len(check.pages)} landscape A4 pages)')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory', type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    combine(args.directory, args.output or args.directory / 'Metatournament_parameter_figures_print.pdf')
