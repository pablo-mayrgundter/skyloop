import { Instances, Instance } from '@react-three/drei'
import { useMemo } from 'react'
import { Vector3 } from 'three'
import { LOS_ANGELES, SAN_FRANCISCO } from '../../config/simulation'
import { buildGreatCircleArc } from '../Globe'

type LoopSectionsProps = {
  count?: number
  elevationKm?: number
}

const LoopSections = ({ count = 12, elevationKm = 20 }: LoopSectionsProps) => {
  const points = useMemo(() => {
    const arc = buildGreatCircleArc(LOS_ANGELES, SAN_FRANCISCO, count)
    const elevated = arc.map((p) => p.clone().add(new Vector3(0, elevationKm * 0.001, 0)))
    return elevated
  }, [count, elevationKm])

  return (
    <Instances limit={points.length} castShadow receiveShadow>
      <sphereGeometry args={[0.01, 16, 16]} />
      <meshStandardMaterial color="#f5c542" emissive="#6b4e00" />
      {points.map((p, idx) => (
        <Instance key={idx} position={p} />
      ))}
    </Instances>
  )
}

export default LoopSections
