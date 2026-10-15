import os
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import image

# === Directory with PNG frames ===
# frame_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "data", "image"))
# output_path = os.path.join(frame_dir, "fluid_bacterium_movie.mp4")

frame_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "data", "image", "velocity_frames"))
frame_files = sorted([f for f in os.listdir(frame_dir) if f.endswith(".png")])
output_path = os.path.join(frame_dir, "tiger_velocity_loop.mp4")

# === Load one frame to get dimensions ===
sample_img = image.imread(os.path.join(frame_dir, frame_files[0]))
fig, ax = plt.subplots(figsize=(10, 8))
img = ax.imshow(sample_img)
ax.axis("off")

# === Frame update function ===
def update(i):
    img.set_data(image.imread(os.path.join(frame_dir, frame_files[i])))
    return [img]

# === Animate ===
ani = animation.FuncAnimation(fig, update, frames=len(frame_files), interval=100)

# === Save as MP4 ===
from matplotlib.animation import FFMpegWriter
writer = FFMpegWriter(fps=10, metadata=dict(artist="AILV"), bitrate=1800)
ani.save(output_path, writer=writer, dpi=200)

plt.close()
print(f"🎬 Movie saved as {output_path}")
