import {GreatCirclePoint} from '../types'

export const EARTH_RADIUS_KM = 6371
export const RENDER_SCALE = 0.001
export const DEFAULT_SIM_SPEED = 1

export const LOS_ANGELES: GreatCirclePoint = { lat: 34.0522, lon: -118.2437 }
export const SAN_FRANCISCO: GreatCirclePoint = { lat: 37.7749, lon: -122.4194 }

export const DEFAULT_ATMOSPHERIC_LAYERS = [
  { name: 'Troposphere', altitudeKm: 12, densityKgPerM3: 1.225 },
  { name: 'Stratosphere', altitudeKm: 50, densityKgPerM3: 0.018 },
  { name: 'Mesosphere', altitudeKm: 85, densityKgPerM3: 0.002 }
]

export const DEFAULT_WIND_PROFILE = [
  { altitudeKm: 0, speedMps: 5, directionDeg: 270 },
  { altitudeKm: 10, speedMps: 12, directionDeg: 280 },
  { altitudeKm: 20, speedMps: 20, directionDeg: 300 }
]
