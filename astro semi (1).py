import os
import time
from tkinter import Tk, filedialog

from astropy.io import fits
from astropy.time import Time

import matplotlib.pyplot as plt
import numpy as np

from openpyxl import Workbook
from openpyxl.styles import Font, Alignment
from openpyxl.utils import get_column_letter


# ============================================================
# ENTER NOAA NUMBER
# ============================================================

noaa_no = input("Enter NOAA Number : ").strip()

if noaa_no == "":
    noaa_no = "Unknown"


# ============================================================
# SELECT FITS FILES
# ============================================================

root = Tk()
root.withdraw()
root.attributes("-topmost", True)

fits_files = filedialog.askopenfilenames(
    title="Select FITS Files",
    filetypes=[
        ("FITS files", "*.fits"),
        ("All files", "*.*")
    ]
)

root.destroy()

fits_files = list(fits_files)

if len(fits_files) < 1:

    print("\nNo FITS files selected. Exiting.")
    input("Press Enter to close...")
    raise SystemExit


print("\nNumber of FITS files selected :", len(fits_files))


# ============================================================
# STORE ALL SELECTED POINTS
# ============================================================

all_points = []


# ============================================================
# PROCESS EACH FITS FILE
# ============================================================

for file_index, fits_file in enumerate(fits_files):

    print("\n" + "=" * 80)
    print(
        f"FILE {file_index + 1} OF {len(fits_files)}"
    )
    print(
        "OPENING:",
        os.path.basename(fits_file)
    )
    print("=" * 80)


    # ========================================================
    # OPEN FITS
    # ========================================================

    try:

        with fits.open(
            fits_file,
            memmap=False
        ) as hdul:

            # Prefer first image HDU
            hdu_idx = 1 if len(hdul) > 1 else 0

            data = hdul[hdu_idx].data

            header = hdul[hdu_idx].header

            if data is None:

                print(
                    "WARNING: No image data found."
                )

                continue

            data = np.asarray(
                data,
                dtype=np.float32
            )


    except Exception as e:

        print(
            "ERROR OPENING FITS FILE:"
        )

        print(e)

        continue


    # ========================================================
    # CHECK IMAGE DIMENSION
    # ========================================================

    if data.ndim != 2:

        print(
            "WARNING: FITS image is not 2-dimensional."
        )

        print(
            "Image shape =",
            data.shape
        )

        continue


    # ========================================================
    # CLEAN DATA
    # ========================================================

    data = np.nan_to_num(
        data,
        nan=0.0,
        posinf=0.0,
        neginf=0.0
    )


    # ========================================================
    # IMPORTANT:
    # DO NOT FLIP HMI DATA
    # ========================================================

    # data = np.fliplr(data)

    # The previous version flipped the image.
    # This can reverse the solar coordinate orientation.
    # We keep the original FITS orientation here.


    # ========================================================
    # READ FITS HEADER
    # ========================================================

    DATE_OBS = header.get(
        "DATE-OBS",
        header.get(
            "DATE_OBS",
            "UNKNOWN"
        )
    )


    CRPIX1 = float(
        header.get(
            "CRPIX1",
            (data.shape[1] + 1) / 2
        )
    )


    CRPIX2 = float(
        header.get(
            "CRPIX2",
            (data.shape[0] + 1) / 2
        )
    )


    CDELT1 = abs(
        float(
            header.get(
                "CDELT1",
                0.5
            )
        )
    )


    RSUN_OBS = float(
        header.get(
            "RSUN_OBS",
            975.0
        )
    )


    # Solar B-angle
    B0 = float(
        header.get(
            "CRLT_OBS",
            header.get(
                "SOLAR_B0",
                0.0
            )
        )
    )


    # P-angle
    P_ANGLE = float(
        header.get(
            "CROTA2",
            0.0
        )
    )


    # ========================================================
    # CONVERT OBSERVATION TIME
    # ========================================================

    try:

        observation_time = Time(
            DATE_OBS
        )


        observation_date = (
            observation_time
            .to_datetime()
            .strftime("%d/%m/%Y")
        )


        observation_datetime = (
            observation_time
            .to_datetime()
            .strftime(
                "%d/%m/%Y %H:%M:%S"
            )
        )


    except Exception:

        observation_date = str(
            DATE_OBS
        )

        observation_datetime = str(
            DATE_OBS
        )


    print(
        "Observation Date :",
        observation_date
    )

    print(
        "Observation Time :",
        observation_datetime
    )


    # ========================================================
    # SOLAR RADIUS IN PIXELS
    # ========================================================

    if CDELT1 <= 0:

        print(
            "WARNING: Invalid CDELT1."
        )

        continue


    solar_radius = (
        RSUN_OBS / CDELT1
    )


    print(
        "Solar radius =",
        f"{solar_radius:.2f}",
        "pixels"
    )


    # ========================================================
    # IMAGE DISPLAY RANGE
    # ========================================================

    finite_data = data[
        np.isfinite(data)
    ]


    if len(finite_data) == 0:

        print(
            "WARNING: Image contains no valid data."
        )

        continue


    low = np.percentile(
        finite_data,
        1
    )


    high = np.percentile(
        finite_data,
        99
    )


    # ========================================================
    # CALCULATE SOLAR COORDINATES
    # ========================================================

    def calculate_values(x, y):

        # ----------------------------------------------------
        # Pixel displacement from disk centre
        # ----------------------------------------------------

        dx = (
            x - CRPIX1
        )

        dy = (
            y - CRPIX2
        )


        r = np.sqrt(
            dx**2 + dy**2
        )


        if r > solar_radius:

            return None


        # ----------------------------------------------------
        # P-ANGLE ROTATION
        # ----------------------------------------------------

        theta = np.radians(
            P_ANGLE
        )


        x_rot = (
            dx * np.cos(theta)
            - dy * np.sin(theta)
        )


        y_rot = (
            dx * np.sin(theta)
            + dy * np.cos(theta)
        )


        # ----------------------------------------------------
        # NORMALIZED SOLAR COORDINATES
        # ----------------------------------------------------

        x_norm = (
            x_rot / solar_radius
        )


        y_norm = (
            y_rot / solar_radius
        )


        r2 = (
            x_norm**2
            + y_norm**2
        )


        if r2 > 1:

            return None


        # ----------------------------------------------------
        # Z COORDINATE
        # ----------------------------------------------------

        z = np.sqrt(
            max(
                0,
                1 - r2
            )
        )


        B0_rad = np.radians(
            B0
        )


        # ----------------------------------------------------
        # HELIOGRAPHIC LATITUDE
        # ----------------------------------------------------

        sin_lat = (
            y_norm * np.cos(B0_rad)
            + z * np.sin(B0_rad)
        )


        sin_lat = np.clip(
            sin_lat,
            -1,
            1
        )


        latitude = np.degrees(
            np.arcsin(
                sin_lat
            )
        )


        # ----------------------------------------------------
        # HELIOGRAPHIC LONGITUDE
        # ----------------------------------------------------

        lon_rad = np.arctan2(
            x_norm,
            (
                z * np.cos(B0_rad)
                - y_norm * np.sin(B0_rad)
            )
        )


        longitude = (
            -np.degrees(
                lon_rad
            )
        )


        # ----------------------------------------------------
        # SOLAR X AND Y
        # ----------------------------------------------------

        solar_x = (
            dx * CDELT1
        )


        solar_y = (
            -dy * CDELT1
        )


        # ----------------------------------------------------
        # NORTH / SOUTH
        # ----------------------------------------------------

        lat_dir = (
            "N"
            if latitude >= 0
            else "S"
        )


        # ----------------------------------------------------
        # EAST / WEST
        # ----------------------------------------------------

        lon_dir = (
            "W"
            if longitude >= 0
            else "E"
        )


        # ----------------------------------------------------
        # COORDINATE TEXT
        # Example: N20W15
        # ----------------------------------------------------

        coord_text = (
            f"{lat_dir}"
            f"{abs(latitude):02.0f}"
            f"{lon_dir}"
            f"{abs(longitude):02.0f}"
        )


        return (
            latitude,
            longitude,
            solar_x,
            solar_y,
            coord_text
        )


    # ========================================================
    # CREATE IMAGE FIGURE
    # ========================================================

    fig, ax = plt.subplots(
        figsize=(9, 9)
    )


    ax.imshow(
        data,
        cmap="gray",
        origin="lower",
        vmin=low,
        vmax=high,
        interpolation="none"
    )


    # ========================================================
    # SOLAR DISK
    # ========================================================

    circle = plt.Circle(
        (
            CRPIX1,
            CRPIX2
        ),
        solar_radius,
        fill=False,
        linewidth=2
    )


    ax.add_patch(
        circle
    )


    # ========================================================
    # ARCSECOND TICKS
    # ========================================================

    ticks_arcsec = np.arange(
        -1000,
        1001,
        500
    )


    tick_pixels_x = (
        ticks_arcsec / CDELT1
    ) + CRPIX1


    tick_pixels_y = (
        -ticks_arcsec / CDELT1
    ) + CRPIX2


    ax.set_xticks(
        tick_pixels_x
    )


    ax.set_yticks(
        tick_pixels_y
    )


    ax.set_xticklabels(
        [
            f"{int(t)}"
            for t in ticks_arcsec
        ]
    )


    ax.set_yticklabels(
        [
            f"{int(t)}"
            for t in ticks_arcsec
        ]
    )


    ax.set_xlabel(
        "Solar X (arcsec)"
    )


    ax.set_ylabel(
        "Solar Y (arcsec)"
    )


    ax.set_title(
        f"{os.path.basename(fits_file)}\n"
        f"Observation Time: "
        f"{observation_datetime}"
    )


    # ========================================================
    # MOUSE COORDINATE DISPLAY
    # ========================================================

    left_display = ax.text(
        0.02,
        0.98,
        "",
        transform=ax.transAxes,
        fontsize=12,
        verticalalignment="top"
    )


    right_display = ax.text(
        0.98,
        0.98,
        "",
        transform=ax.transAxes,
        fontsize=12,
        verticalalignment="top",
        horizontalalignment="right"
    )


    last_update = [0]


    # ========================================================
    # MOUSE MOVE
    # ========================================================

    def mouse_move(event):

        current_time = time.time()


        if (
            current_time
            - last_update[0]
            < 0.03
        ):

            return


        last_update[0] = current_time


        if (
            event.xdata is None
            or event.ydata is None
        ):

            return


        result = calculate_values(
            event.xdata,
            event.ydata
        )


        if result is None:

            left_display.set_text("")
            right_display.set_text("")

            fig.canvas.draw_idle()

            return


        (
            latitude,
            longitude,
            solar_x,
            solar_y,
            coord_text
        ) = result


        left_display.set_text(
            f"{coord_text}\n"
            f'X = {solar_x:.0f}"\n'
            f'Y = {solar_y:.0f}"'
        )


        right_display.set_text(
            f"Latitude = "
            f"{latitude:.2f}°\n"
            f"Longitude = "
            f"{longitude:.2f}°"
        )


        fig.canvas.draw_idle()


    # ========================================================
    # MOUSE CLICK
    # ========================================================

    def onclick(event):

        if (
            event.xdata is None
            or event.ydata is None
        ):

            return


        result = calculate_values(
            event.xdata,
            event.ydata
        )


        if result is None:

            return


        (
            latitude,
            longitude,
            solar_x,
            solar_y,
            coord_text
        ) = result


        # ----------------------------------------------------
        # STORE POINT
        # ----------------------------------------------------

        current_point = {

            "file":
                os.path.basename(
                    fits_file
                ),

            "time":
                DATE_OBS,

            "date":
                observation_date,

            "datetime":
                observation_datetime,

            "latitude":
                latitude,

            "longitude":
                longitude,

            "solar_x":
                solar_x,

            "solar_y":
                solar_y,

            "coord":
                coord_text
        }


        all_points.append(
            current_point
        )


        # ----------------------------------------------------
        # PRINT SELECTED POINT
        # ----------------------------------------------------

        print(
            "\nSaved Point -> "
            f"Date: {observation_date} | "
            f"Coord: {coord_text} | "
            f'X: {solar_x:.0f}" | '
            f'Y: {solar_y:.0f}" | '
            f"Lat: {latitude:.2f}° | "
            f"Lon: {longitude:.2f}°"
        )


        # ----------------------------------------------------
        # SHOW POINT ON IMAGE
        # ----------------------------------------------------

        ax.scatter(
            event.xdata,
            event.ydata,
            s=90,
            facecolors="red",
            edgecolors="white",
            linewidths=1.5,
            alpha=0.85
        )


        ax.text(
            event.xdata + 20,
            event.ydata + 20,
            coord_text,
            fontsize=10
        )


        fig.canvas.draw_idle()


    # ========================================================
    # CONNECT MOUSE EVENTS
    # ========================================================

    fig.canvas.mpl_connect(
        "motion_notify_event",
        mouse_move
    )


    fig.canvas.mpl_connect(
        "button_press_event",
        onclick
    )


    # ========================================================
    # SHOW IMAGE
    # ========================================================

    plt.show()


    # ========================================================
    # CLOSE FIGURE AFTER SELECTION
    # ========================================================

    plt.close(
        fig
    )


# ============================================================
# CHECK NUMBER OF SELECTED POINTS
# ============================================================

print("\n" + "=" * 80)

print(
    "TOTAL SELECTED POINTS =",
    len(all_points)
)

print("=" * 80)


if len(all_points) < 2:

    print(
        "\nNot enough points selected."
    )

    print(
        "Minimum 2 points are required "
        "to calculate rotation."
    )

    input(
        "\nPress Enter to close..."
    )

    raise SystemExit


# ============================================================
# SORT POINTS BY OBSERVATION TIME
# ============================================================

try:

    all_points.sort(
        key=lambda p:
        Time(
            p["time"]
        ).jd
    )

except Exception as e:

    print(
        "ERROR SORTING OBSERVATION TIMES:"
    )

    print(e)

    input(
        "\nPress Enter to close..."
    )

    raise SystemExit


# ============================================================
# CALCULATE ROTATION
# ============================================================

rotation_rates = []

latitudes_plot = []

rotation_plot = []

rotation_rows = []


for i in range(
    1,
    len(all_points)
):

    p1 = all_points[
        i - 1
    ]

    p2 = all_points[
        i
    ]


    # --------------------------------------------------------
    # TIME
    # --------------------------------------------------------

    try:

        t1 = Time(
            p1["time"]
        )

        t2 = Time(
            p2["time"]
        )


        dt_days = (
            t2 - t1
        ).to_value(
            "day"
        )


    except Exception:

        print(
            "WARNING: Could not calculate time "
            f"difference for points {i} and {i + 1}."
        )

        continue


    if dt_days <= 0:

        print(
            "WARNING: Non-positive time interval."
        )

        continue


    # --------------------------------------------------------
    # LONGITUDE DIFFERENCE
    # --------------------------------------------------------

    dlon = (
        p2["longitude"]
        - p1["longitude"]
    )


    # Handle 360-degree wrap
    if dlon > 180:

        dlon -= 360


    elif dlon < -180:

        dlon += 360


    # --------------------------------------------------------
    # SYNODIC ROTATION RATE
    # --------------------------------------------------------

    synodic_rate = (
        dlon / dt_days
    )


    # --------------------------------------------------------
    # ABSOLUTE RATE FOR ORIGINAL PLOT
    # --------------------------------------------------------

    absolute_synodic_rate = abs(
        synodic_rate
    )


    # --------------------------------------------------------
    # SIDEREAL ROTATION RATE
    #
    # Earth orbital angular motion:
    # 0.9856 degree/day
    #
    # For westward solar rotation:
    #
    # Sidereal = Synodic + 0.9856
    # --------------------------------------------------------

    sidereal_rate = (
        synodic_rate
        + 0.9856
    )


    absolute_sidereal_rate = abs(
        sidereal_rate
    )


    # --------------------------------------------------------
    # AVERAGE LATITUDE
    # --------------------------------------------------------

    avg_latitude = (
        p1["latitude"]
        + p2["latitude"]
    ) / 2


    # --------------------------------------------------------
    # STORE
    # --------------------------------------------------------

    rotation_rates.append(
        absolute_synodic_rate
    )


    latitudes_plot.append(
        avg_latitude
    )


    rotation_plot.append(
        absolute_synodic_rate
    )


    rotation_rows.append(
        {
            "point_number":
                i + 1,

            "from_date":
                p1["datetime"],

            "to_date":
                p2["datetime"],

            "latitude":
                avg_latitude,

            "longitude_change":
                dlon,

            "time_difference":
                dt_days,

            "synodic_rate":
                absolute_synodic_rate,

            "sidereal_rate":
                absolute_sidereal_rate
        }
    )


# ============================================================
# CHECK ROTATION DATA
# ============================================================

if len(rotation_rates) == 0:

    print(
        "\nNo valid rotation rates calculated."
    )

    input(
        "\nPress Enter to close..."
    )

    raise SystemExit


# ============================================================
# AVERAGE SYNODIC ROTATION RATE
# ============================================================

avg_synodic_rotation = np.mean(
    rotation_rates
)


# ============================================================
# AVERAGE SIDEREAL ROTATION RATE
# ============================================================

earth_orbital_rate = 0.9856


avg_sidereal_rotation = (
    avg_synodic_rotation
    + earth_orbital_rate
)


# ============================================================
# ROTATION PERIODS
# ============================================================

if avg_synodic_rotation > 0:

    synodic_period = (
        360
        / avg_synodic_rotation
    )

else:

    synodic_period = np.nan


if avg_sidereal_rotation > 0:

    sidereal_period = (
        360
        / avg_sidereal_rotation
    )

else:

    sidereal_period = np.nan


# ============================================================
# FINAL ROTATION OUTPUT
# ============================================================

print("\n" + "=" * 80)

print(
    "AVERAGE SOLAR ROTATION RESULTS"
)

print("=" * 80)

print(
    f"Average Synodic Rotation Rate  = "
    f"{avg_synodic_rotation:.4f} deg/day"
)

print(
    f"Average Sidereal Rotation Rate = "
    f"{avg_sidereal_rotation:.4f} deg/day"
)

print()

print(
    f"Average Synodic Rotation Period  = "
    f"{synodic_period:.4f} days"
)

print(
    f"Average Sidereal Rotation Period = "
    f"{sidereal_period:.4f} days"
)

print("=" * 80)


# ============================================================
# DIFFERENTIAL ROTATION GRAPH
# ============================================================

fig4, ax4 = plt.subplots(
    figsize=(8, 6)
)


# ------------------------------------------------------------
# SORT BY LATITUDE
# ------------------------------------------------------------

plot_data = sorted(
    zip(
        latitudes_plot,
        rotation_plot
    )
)


sorted_latitudes = [
    item[0]
    for item in plot_data
]


sorted_rotation = [
    item[1]
    for item in plot_data
]


# ------------------------------------------------------------
# PLOT
# ------------------------------------------------------------

ax4.plot(
    sorted_latitudes,
    sorted_rotation,
    marker="o",
    linestyle="-",
    linewidth=1.5,
    markersize=6
)


ax4.set_title(
    "Solar Differential Rotation"
)


ax4.set_xlabel(
    "Latitude (degree)"
)


ax4.set_ylabel(
    "Synodic Rotation Rate (degree/day)"
)


ax4.grid(
    True,
    linestyle="--",
    alpha=0.5
)


fig4.tight_layout()


# ============================================================
# SAVE GRAPH
# ============================================================

plot_filename = (
    f"NOAA_{noaa_no}_"
    f"solar_rotation_plot.png"
)


try:

    fig4.savefig(
        plot_filename,
        dpi=300,
        bbox_inches="tight"
    )

    print(
        "\nGraph saved successfully:"
    )

    print(
        os.path.abspath(
            plot_filename
        )
    )


except Exception as e:

    print(
        "\nERROR SAVING GRAPH:"
    )

    print(e)


plt.show()

plt.close(
    fig4
)


# ============================================================
# CREATE EXCEL WORKBOOK
# ============================================================

try:

    wb = Workbook()


    ws = wb.active


    ws.title = (
        "Sunspot Rotation Data"
    )


    # ========================================================
    # EXCEL HEADERS
    # ========================================================

    headers = [

        "NOAA",

        "DATE",

        "DATE & TIME",

        "FITS FILE",

        "COORDINATES",

        "LATITUDE (deg)",

        "LONGITUDE (deg)",

        "SOLAR X (arcsec)",

        "SOLAR Y (arcsec)",

        "SYNODIC ROTATION RATE (deg/day)",

        "SIDEREAL ROTATION RATE (deg/day)"
    ]


    # --------------------------------------------------------
    # WRITE HEADERS
    # --------------------------------------------------------

    for col_num, header_name in enumerate(
        headers,
        start=1
    ):

        cell = ws.cell(
            row=1,
            column=col_num,
            value=header_name
        )


        cell.font = Font(
            bold=True
        )


        cell.alignment = Alignment(
            horizontal="center",
            vertical="center"
        )


    # ========================================================
    # WRITE POINT DATA
    # ========================================================

    for i, point in enumerate(
        all_points
    ):

        row_num = i + 2


        # ----------------------------------------------------
        # NOAA
        # ----------------------------------------------------

        ws.cell(
            row=row_num,
            column=1,
            value=noaa_no
        )


        # ----------------------------------------------------
        # DATE AND TIME
        # ----------------------------------------------------

        try:

            excel_date = (
                Time(
                    point["time"]
                )
                .to_datetime()
            )


            date_cell = ws.cell(
                row=row_num,
                column=2,
                value=excel_date
            )


            date_cell.number_format = (
                "dd/mm/yyyy"
            )


            datetime_cell = ws.cell(
                row=row_num,
                column=3,
                value=excel_date
            )


            datetime_cell.number_format = (
                "dd/mm/yyyy hh:mm:ss"
            )


        except Exception:

            ws.cell(
                row=row_num,
                column=2,
                value=point["date"]
            )


            ws.cell(
                row=row_num,
                column=3,
                value=point["datetime"]
            )


        # ----------------------------------------------------
        # FITS FILE
        # ----------------------------------------------------

        ws.cell(
            row=row_num,
            column=4,
            value=point["file"]
        )


        # ----------------------------------------------------
        # COORDINATES
        # ----------------------------------------------------

        ws.cell(
            row=row_num,
            column=5,
            value=point["coord"]
        )


        # ----------------------------------------------------
        # LATITUDE
        # ----------------------------------------------------

        ws.cell(
            row=row_num,
            column=6,
            value=round(
                point["latitude"],
                4
            )
        )


        # ----------------------------------------------------
        # LONGITUDE
        # ----------------------------------------------------

        ws.cell(
            row=row_num,
            column=7,
            value=round(
                point["longitude"],
                4
            )
        )


        # ----------------------------------------------------
        # SOLAR X
        # ----------------------------------------------------

        ws.cell(
            row=row_num,
            column=8,
            value=round(
                point["solar_x"],
                2
            )
        )


        # ----------------------------------------------------
        # SOLAR Y
        # ----------------------------------------------------

        ws.cell(
            row=row_num,
            column=9,
            value=round(
                point["solar_y"],
                2
            )
        )


        # ----------------------------------------------------
        # ROTATION RATES
        # ----------------------------------------------------

        if i == 0:

            ws.cell(
                row=row_num,
                column=10,
                value="-"
            )


            ws.cell(
                row=row_num,
                column=11,
                value="-"
            )


        else:

            p1 = all_points[
                i - 1
            ]


            p2 = all_points[
                i
            ]


            try:

                t1 = Time(
                    p1["time"]
                )


                t2 = Time(
                    p2["time"]
                )


                dt_days = (
                    t2 - t1
                ).to_value(
                    "day"
                )


            except Exception:

                dt_days = 0


            if dt_days > 0:

                dlon = (
                    p2["longitude"]
                    - p1["longitude"]
                )


                # 360 degree correction
                if dlon > 180:

                    dlon -= 360


                elif dlon < -180:

                    dlon += 360


                # --------------------------------------------
                # SYNODIC
                # --------------------------------------------

                synodic_rate = abs(
                    dlon / dt_days
                )


                # --------------------------------------------
                # SIDEREAL
                # --------------------------------------------

                sidereal_rate = (
                    synodic_rate
                    + earth_orbital_rate
                )


                ws.cell(
                    row=row_num,
                    column=10,
                    value=round(
                        synodic_rate,
                        4
                    )
                )


                ws.cell(
                    row=row_num,
                    column=11,
                    value=round(
                        sidereal_rate,
                        4
                    )
                )


            else:

                ws.cell(
                    row=row_num,
                    column=10,
                    value="-"
                )


                ws.cell(
                    row=row_num,
                    column=11,
                    value="-"
                )


    # ========================================================
    # SUMMARY
    # ========================================================

    summary_row = (
        len(all_points)
        + 4
    )


    # --------------------------------------------------------
    # SUMMARY TITLE
    # --------------------------------------------------------

    ws.cell(
        row=summary_row,
        column=1,
        value="SUMMARY"
    )


    ws.cell(
        row=summary_row,
        column=1
    ).font = Font(
        bold=True
    )


    # --------------------------------------------------------
    # EARTH ORBITAL RATE
    # --------------------------------------------------------

    ws.cell(
        row=summary_row + 1,
        column=1,
        value="Earth Orbital Rate"
    )


    ws.cell(
        row=summary_row + 1,
        column=2,
        value=0.9856
    )


    ws.cell(
        row=summary_row + 1,
        column=3,
        value="deg/day"
    )


    # --------------------------------------------------------
    # AVERAGE SYNODIC RATE
    # --------------------------------------------------------

    ws.cell(
        row=summary_row + 2,
        column=1,
        value="Average Synodic Rotation Rate"
    )


    ws.cell(
        row=summary_row + 2,
        column=2,
        value=round(
            avg_synodic_rotation,
            4
        )
    )


    ws.cell(
        row=summary_row + 2,
        column=3,
        value="deg/day"
    )


    # --------------------------------------------------------
    # AVERAGE SIDEREAL RATE
    # --------------------------------------------------------

    ws.cell(
        row=summary_row + 3,
        column=1,
        value="Average Sidereal Rotation Rate"
    )


    ws.cell(
        row=summary_row + 3,
        column=2,
        value=round(
            avg_sidereal_rotation,
            4
        )
    )


    ws.cell(
        row=summary_row + 3,
        column=3,
        value="deg/day"
    )


    # --------------------------------------------------------
    # SYNODIC PERIOD
    # --------------------------------------------------------

    ws.cell(
        row=summary_row + 4,
        column=1,
        value="Average Synodic Rotation Period"
    )


    ws.cell(
        row=summary_row + 4,
        column=2,
        value=round(
            synodic_period,
            4
        )
    )


    ws.cell(
        row=summary_row + 4,
        column=3,
        value="days"
    )


    # --------------------------------------------------------
    # SIDEREAL PERIOD
    # --------------------------------------------------------

    ws.cell(
        row=summary_row + 5,
        column=1,
        value="Average Sidereal Rotation Period"
    )


    ws.cell(
        row=summary_row + 5,
        column=2,
        value=round(
            sidereal_period,
            4
        )
    )


    ws.cell(
        row=summary_row + 5,
        column=3,
        value="days"
    )


    # ========================================================
    # COLUMN WIDTH
    # ========================================================

    column_widths = {

        1: 12,

        2: 15,

        3: 22,

        4: 45,

        5: 18,

        6: 18,

        7: 19,

        8: 20,

        9: 20,

        10: 35,

        11: 35
    }


    for col_num, width in (
        column_widths.items()
    ):

        ws.column_dimensions[
            get_column_letter(
                col_num
            )
        ].width = width


    # ========================================================
    # ALIGNMENT
    # ========================================================

    for row in ws.iter_rows():

        for cell in row:

            cell.alignment = Alignment(
                vertical="center",
                horizontal="center"
            )


    # ========================================================
    # FREEZE HEADER
    # ========================================================

    ws.freeze_panes = "A2"


    # ========================================================
    # AUTO FILTER
    # ========================================================

    ws.auto_filter.ref = (
        f"A1:K{len(all_points) + 1}"
    )


    # ========================================================
    # SAVE EXCEL
    # ========================================================

    output_excel = (
        f"NOAA_{noaa_no}_"
        f"Sunspot_Rotation.xlsx"
    )


    # --------------------------------------------------------
    # USE ABSOLUTE PATH
    # --------------------------------------------------------

    output_excel = os.path.abspath(
        output_excel
    )


    wb.save(
        output_excel
    )


    # ========================================================
    # VERIFY EXCEL FILE
    # ========================================================

    if os.path.exists(
        output_excel
    ):

        print(
            "\nExcel file saved successfully:"
        )

        print(
            output_excel
        )

        print(
            "File size:",
            os.path.getsize(
                output_excel
            ),
            "bytes"
        )


    else:

        print(
            "\nWARNING:"
        )

        print(
            "Excel file was not found after saving."
        )


except Exception as e:

    print(
        "\n" + "=" * 80
    )

    print(
        "ERROR CREATING EXCEL FILE"
    )

    print(
        "=" * 80
    )

    print(
        repr(e)
    )

    print(
        "\nMake sure the Excel file is not already "
        "open in Microsoft Excel."
    )


# ============================================================
# FINAL OUTPUT
# ============================================================

print("\n" + "=" * 80)

print(
    "FINAL RESULTS"
)

print("=" * 80)

print(
    f"NOAA Number = {noaa_no}"
)

print(
    f"Number of selected points = "
    f"{len(all_points)}"
)

print(
    f"Number of rotation intervals = "
    f"{len(rotation_rates)}"
)

print()

print(
    f"Average Synodic Rotation Rate = "
    f"{avg_synodic_rotation:.4f} deg/day"
)

print(
    f"Average Sidereal Rotation Rate = "
    f"{avg_sidereal_rotation:.4f} deg/day"
)

print()

print(
    f"Average Synodic Rotation Period = "
    f"{synodic_period:.4f} days"
)

print(
    f"Average Sidereal Rotation Period = "
    f"{sidereal_period:.4f} days"
)

print()

print(
    "Graph file =",
    os.path.abspath(
        plot_filename
    )
)

print(
    "Excel file =",
    os.path.abspath(
        output_excel
    )
)

print("=" * 80)

print(
    "\nDONE"
)

input(
    "\nPress Enter to close..."
)
