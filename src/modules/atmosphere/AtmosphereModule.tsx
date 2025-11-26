import { AtmosphericLayer } from '../../types/simulation'

const layers: AtmosphericLayer[] = [
  { name: 'Troposphere', altitudeKm: 12, density: 1.225 },
  { name: 'Stratosphere', altitudeKm: 50, density: 0.1 },
  { name: 'Mesosphere', altitudeKm: 85, density: 0.003 }
]

const AtmosphereModule = () => {
  return (
    <div className="controls-grid">
      <h3>Atmospheric layers</h3>
      {layers.map((layer) => (
        <div key={layer.name} className="metric">
          <strong>{layer.name}</strong>
          <div>Altitude: {layer.altitudeKm} km</div>
          <div>Density: {layer.density} kg/m³</div>
        </div>
      ))}
      <p>TODO: Integrate atmospheric drag, heating, and refraction models.</p>
    </div>
  )
}

export default AtmosphereModule
