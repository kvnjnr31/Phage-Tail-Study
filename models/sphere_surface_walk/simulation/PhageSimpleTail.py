import numpy as np

class PhageSimpleTail:
    """
    Minimal bacteriophage model using a head and a simple linear tail (beads-on-a-stick).

    Attributes:
    - head_pos: 3D position of the phage head
    - head_radius: size of the phage head (sphere)
    - tail_length: number of segments (beads) in the tail
    - tail_spacing: distance between each bead
    - tail_coords: 3D coordinates of each bead (tail)
    - bound: boolean indicating whether phage is bound to a bacterium
    """

    def __init__(self, head_pos=np.array([30.0, 15.0, 15.0]), head_radius=1.5,
                 tail_length=10, tail_spacing=1.0):
        self.head_pos = np.array(head_pos, dtype=float)
        self.head_radius = head_radius
        self.tail_length = tail_length
        self.tail_spacing = tail_spacing
        self.tail_coords = self.initialize_tail()
        self.bound = False

    def initialize_tail(self):
        """Lay out the tail beads in a straight line extending from the head."""
        tail = []
        for i in range(1, self.tail_length + 1):
            bead = self.head_pos - np.array([i * self.tail_spacing, 0.0, 0.0])
            tail.append(bead)
        return np.array(tail)

    def update(self, flow_field, dt=0.1):
        """
        Optional dynamics step (placeholder): advect head or allow tail to wobble.
        In minimal version, phage is stationary or moves rigidly.
        """
        pass

    def check_binding(self, bacterium_pos, threshold=1.0):
        """
        Check if any tail bead is close enough to bind to the bacterium.
        """
        for bead in self.tail_coords:
            if np.linalg.norm(bead - bacterium_pos) < threshold:
                self.bound = True
                return True
        return False
