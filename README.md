# pointdexter

It's just connect the dots. No, it's not impossible.


### Inputs:

#### Detection Catalog 

Provide a CSV file containing the detections to be linked. The CSV file should have the following columns:

| Column Name   | Description                         | Units        | Required      |
|---------------|-------------------------------------|--------------|---------------|
| `ra`          | Right Ascension (Equatorial)        | degrees      | Yes           |
| `dec`         | Declination (Equatorial)            | degrees      | Yes           |
| `epoch`       | Epoch of observation                | UTC JD       | Yes           |
| `obscode`     | Observatory code                    | string       | If no obs pos |
| `obs_x`       | Observer X position (Equatorial)    | au           | If no obscode |
| `obs_y`       | Observer Y position (Equatorial)    | au           | If no obscode |
| `obs_z`       | Observer Z position (Equatorial)    | au           | If no obscode |
| `obs_vx`      | Observer X velocity (Equatorial)    | au/day       | No            |
| `obs_vy`      | Observer Y velocity (Equatorial)    | au/day       | No            |
| `obs_vz`      | Observer Z velocity (Equatorial)    | au/day       | No            |
| `ast_ucty`    | Astrometric uncertainty             | arcseconds   | No            |
| `magnitude`   | Apparent magnitude of detection     | magnitudes   | No            |
| `filter`      | Photometric filter                  | string       | No            |
| `detid`       | Detection ID                        | string       | No            |
| `trackid`     | Tracklet ID                         | string       | No            |

#### Initial Conditions File

Provide a CSV file containing the initial conditions for orbit propagation. There are two supported formats: "spherical" and "keplerian". Specify the desired format using the `ic_type` parameter in the configuration file.

The config file should contain the keyword `reference_epoch`, which specifies the epoch (in UTC JD) at which the initial conditions are defined. 
It should also contain the keyword `ic_origin`, which specifies the origin of the initial conditions, either "ssb" (Solar System Barycenter) or "sun".

##### Spherical Initial Conditions
| Column Name   | Description                         | Units        | Required      |
|---------------|-------------------------------------|--------------|---------------|
| `r`           | Heliocentric distance               | au           | Yes           |
| `vr`         | Radial velocity                    | au/day       | Yes           |
| `vo`     | Transverse velocity                 | au/day       | Yes           |
| `inc`         | Inclination                         | radians      | Yes           |
| `kappa`       | The parity of the inclination        | unitless $\pm 1$     | Yes           |
| `id`          | Initial condition ID                | string       | Yes            |


##### Keplerian Initial Conditions
| Column Name   | Description                         | Units        | Required      |
|---------------|-------------------------------------|--------------|---------------|
| `q`           | Perihelion distance                 | au           | Yes           |
| `e`           | Eccentricity                        | unitless     | Yes           |
| `inc`         | Inclination                         | radians      | Yes           |
| `M0`       | Mean anomaly at the reference epoch               | radians      | Yes           |
| `kappa`       | The parity of the inclination        | unitless $\pm 1$     | Yes           |
| `id`          | Initial condition ID                | string       | Yes            |