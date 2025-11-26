import { Vector3 } from 'three'

// Earth radius in km (Three.js uses 1 unit = 1 km)
export const EARTH_RADIUS_KM = 6371 // 6,371 km

// Loop altitude above Earth surface in km (~20km altitude)
export const LOOP_ALTITUDE_KM = 20

// Loop radius from Earth center in km
export const LOOP_RADIUS_KM = EARTH_RADIUS_KM + LOOP_ALTITUDE_KM

export const LOS_ANGELES = {
  lat: 34.0522,
  lon: -118.2437
}

export const SAN_FRANCISCO = {
  lat: 37.7749,
  lon: -122.4194
}

export const DEFAULT_SIM_SPEED = 1

export const greatCirclePoint = (lat: number, lon: number, radius = EARTH_RADIUS_KM) => {
  const phi = (90 - lat) * (Math.PI / 180)
  const theta = (lon + 180) * (Math.PI / 180)
  const x = -radius * Math.sin(phi) * Math.cos(theta)
  const z = radius * Math.sin(phi) * Math.sin(theta)
  const y = radius * Math.cos(phi)
  return new Vector3(x, y, z)
}

export const interpolateGreatCircle = (
  start: { lat: number; lon: number },
  end: { lat: number; lon: number },
  t: number,
  radius = EARTH_RADIUS_KM
) => {
  const startVec = greatCirclePoint(start.lat, start.lon, radius).normalize()
  const endVec = greatCirclePoint(end.lat, end.lon, radius).normalize()
  const omega = Math.acos(Math.max(-1, Math.min(1, startVec.dot(endVec))))
  if (omega === 0) {
    return startVec.clone().multiplyScalar(radius)
  }
  const sinOmega = Math.sin(omega)
  const factorA = Math.sin((1 - t) * omega) / sinOmega
  const factorB = Math.sin(t * omega) / sinOmega
  return startVec.clone().multiplyScalar(factorA * radius).add(endVec.multiplyScalar(factorB * radius))
}
