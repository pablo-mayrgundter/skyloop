import { DEFAULT_ATMOSPHERIC_LAYERS } from '../../config/simulation'

const AtmosphereModule = () => (
  <div className="panel">
    <h3>Atmosphere Lab</h3>
    <p>Profile atmospheric density and drag interactions.</p>
    <ul>
      {DEFAULT_ATMOSPHERIC_LAYERS.map((layer) => (
        <li key={layer.name}>
          <strong>{layer.name}</strong>: {layer.altitudeKm} km, {layer.densityKgPerM3} kg/m³
        </li>
      ))}
    </ul>
    <p>TODO: integrate rayleigh scattering shader + drag estimator.</p>
  </div>
)

export default AtmosphereModule
