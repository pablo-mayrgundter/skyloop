import { DEFAULT_WIND_PROFILE } from '../../config/simulation'

const WeatherModule = () => (
  <div className="panel">
    <h3>Weather Lab</h3>
    <p>Wind vectors drive lateral load cases.</p>
    <table>
      <thead>
        <tr>
          <th>Altitude (km)</th>
          <th>Speed (m/s)</th>
          <th>Direction (deg)</th>
        </tr>
      </thead>
      <tbody>
        {DEFAULT_WIND_PROFILE.map((wind) => (
          <tr key={wind.altitudeKm}>
            <td>{wind.altitudeKm}</td>
            <td>{wind.speedMps}</td>
            <td>{wind.directionDeg}</td>
          </tr>
        ))}
      </tbody>
    </table>
    <p>TODO: connect to NOAA datasets or procedural wind fields.</p>
  </div>
)

export default WeatherModule
