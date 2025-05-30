from enum import Enum, unique
from typing import Optional
from sympy import acos, atan2, cos, sin, sqrt, Expr
from sympy.vector import CoordSys3D, Vector as SymVector
from ..symbols.symbols import next_name


class CoordinateSystem:
    __slots__ = "_coord_system", "_coord_system_type"

    _coord_system: CoordSys3D
    _coord_system_type: "System"

    class NotSupported(ValueError):
        """This coordinate system is not supported"""

    @unique
    class System(Enum):
        CARTESIAN = 0
        CYLINDRICAL = 1
        SPHERICAL = 2

        @staticmethod
        def get_name(coord_system_type):
            match coord_system_type:
                case CoordinateSystem.System.CARTESIAN:
                    return "cartesian"
                case CoordinateSystem.System.CYLINDRICAL:
                    return "cylindrical"
                case CoordinateSystem.System.SPHERICAL:
                    return "spherical"
                case _:
                    raise CoordinateSystem.NotSupported(coord_system_type)

    @staticmethod
    def system_to_transformation_name(coord_system_type: System) -> str:
        return CoordinateSystem.System.get_name(coord_system_type)

    @staticmethod
    def system_to_base_scalars(coord_system_type: System) -> tuple[str, str, str]:
        match coord_system_type:
            case CoordinateSystem.System.CARTESIAN:
                return "x", "y", "z"
            case CoordinateSystem.System.CYLINDRICAL:
                return "r", "theta", "phi"
            case CoordinateSystem.System.SPHERICAL:
                # theta - azimuthal angle
                # phi - polar angle
                return "r", "theta", "phi"
            case _:
                raise CoordinateSystem.NotSupported(coord_system_type)

    @staticmethod
    def is_angle_component(coord_system_type: System, component_idx: int) -> bool:
        if coord_system_type == CoordinateSystem.System.CARTESIAN:
            return False
        if coord_system_type == CoordinateSystem.System.CYLINDRICAL:
            return component_idx == 1
        if coord_system_type == CoordinateSystem.System.SPHERICAL:
            return component_idx in (1, 2)
        return False

    def __init__(self,
        coord_system_type: System = System.CARTESIAN,
        inner: Optional[CoordSys3D] = None) -> None:
        self._coord_system_type = coord_system_type
        self._coord_system = inner or CoordSys3D(
            next_name("SYS"),
            variable_names=CoordinateSystem.system_to_base_scalars(coord_system_type)
        )

    @property
    def coord_system(self) -> CoordSys3D:
        return self._coord_system

    @property
    def coord_system_type(self) -> System:
        return self._coord_system_type

    def transformation_to_system(self, coord_system_type: System) -> tuple[Expr, Expr, Expr]:
        if self._coord_system_type == self.System.CYLINDRICAL:
            r, theta, z = self._coord_system.base_scalars()
            cylindrical_conversions = {
                self.System.CARTESIAN: (r * cos(theta), r * sin(theta), z),
                self.System.CYLINDRICAL: (r, theta, z)
            }
            transformation = cylindrical_conversions.get(coord_system_type)
            if transformation is not None:
                return transformation

        if self._coord_system_type == self.System.SPHERICAL:
            r, theta, phi = self._coord_system.base_scalars()
            spherical_conversions = {
                self.System.CARTESIAN:
                (r * cos(theta) * sin(phi), r * sin(theta) * sin(phi), r * cos(phi)),
                self.System.SPHERICAL: (r, theta, phi)
            }
            transformation = spherical_conversions.get(coord_system_type)
            if transformation is not None:
                return transformation

        if self._coord_system_type == self.System.CARTESIAN:
            x, y, z = self._coord_system.base_scalars()
            cartesian_conversions = {
                self.System.CYLINDRICAL: (sqrt(x**2 + y**2), atan2(y, x), z),
                self.System.SPHERICAL: (sqrt(x**2 + y**2 + z**2), atan2(y,
                x), acos(z / sqrt(x**2 + y**2 + z**2))),
                self.System.CARTESIAN: (x, y, z)
            }
            transformation = cartesian_conversions.get(coord_system_type)
            if transformation is not None:
                return transformation

        coord_name_from = self.system_to_transformation_name(self._coord_system_type)
        coord_name_to = self.system_to_transformation_name(coord_system_type)
        raise ValueError(
            f"Transformation is not supported: from {coord_name_from} to {coord_name_to}"
        )


def coordinates_transform(
    from_system: CoordinateSystem,
    coord_system_type: CoordinateSystem.System = CoordinateSystem.System.CARTESIAN,
) -> CoordinateSystem:
    """Change coordinate system type, e.g. from cartesian to cylindrical"""

    new_coord_system = from_system.coord_system.create_new(
        next_name("SYS"),
        variable_names=CoordinateSystem.system_to_base_scalars(coord_system_type),
        transformation=None
    )
    return CoordinateSystem(coord_system_type, new_coord_system)


def coordinates_rotate(self: CoordinateSystem, angle: Expr, axis: SymVector) -> CoordinateSystem:
    if self.coord_system_type != CoordinateSystem.System.CARTESIAN:
        coord_name_from = CoordinateSystem.system_to_transformation_name(self.coord_system_type)
        raise ValueError(f"Rotation only supported for cartesian coordinates, got {coord_name_from}")
    return CoordinateSystem(
        self.coord_system_type,
        self.coord_system.orient_new_axis(next_name("C"), angle, axis)
    )
