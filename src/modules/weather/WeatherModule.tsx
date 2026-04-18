import { WindProfile } from '../../types/simulation'

const WeatherModule = () => {
  const winds: WindProfile[] = [
    { altitudeKm: 0, speed: 5, directionDeg: 260 },
    { altitudeKm: 10, speed: 40, directionDeg: 240 },
    { altitudeKm: 20, speed: 20, directionDeg: 220 }
  ]

  return (
    <div className="controls-grid">
      <h3>Wind profiles</h3>
      {winds.map((wind) => (
        <div key={wind.altitudeKm} className="metric">
          <div>{wind.altitudeKm} km</div>
          <strong>{wind.speed} m/s @ {wind.directionDeg}°</strong>
        </div>
      ))}
      <p>TODO: couple wind profiles into structural loads and guidance.</p>
    </div>
  )
}

export default WeatherModule
