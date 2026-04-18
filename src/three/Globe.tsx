import { Line } from '@react-three/drei'
import { useMemo } from 'react'
import { Vector3 } from 'three'
import { EARTH_RADIUS_KM, LOS_ANGELES, RENDER_SCALE, SAN_FRANCISCO } from '../config/simulation'
import { GreatCirclePoint } from '../types'

const degToRad = (deg: number) => (deg * Math.PI) / 180

export const latLonToVector3 = (point: GreatCirclePoint, radiusKm = EARTH_RADIUS_KM) => {
  const phi = degToRad(90 - point.lat)
  const theta = degToRad(point.lon + 180)
  const r = radiusKm * RENDER_SCALE
  const x = r * Math.sin(phi) * Math.cos(theta)
  const y = r * Math.cos(phi)
  const z = r * Math.sin(phi) * Math.sin(theta)
  return new Vector3(x, y, z)
}

const buildGreatCircleArc = (start: GreatCirclePoint, end: GreatCirclePoint, segments = 64) => {
  const startVec = latLonToVector3(start)
  const endVec = latLonToVector3(end)
  const points: Vector3[] = []
  for (let i = 0; i <= segments; i += 1) {
    const t = i / segments
    const interpolated = new Vector3().lerpVectors(startVec, endVec, t).normalize().multiplyScalar(startVec.length())
    points.push(interpolated)
  }
  return points
}

const Globe = () => {
  const arcPoints = useMemo(() => buildGreatCircleArc(LOS_ANGELES, SAN_FRANCISCO), [])
  const radius = EARTH_RADIUS_KM * RENDER_SCALE
  return (
    <group>
      <mesh>
        <sphereGeometry args={[radius, 64, 64]} />
        <meshStandardMaterial color="#1b3a5a" emissive="#0a1a30" wireframe={false} />
      </mesh>
      <Line points={arcPoints} color="#78e0ff" lineWidth={2} />
    </group>
  )
}

export default Globe
export { buildGreatCircleArc }
