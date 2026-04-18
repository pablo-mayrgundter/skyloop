import { useMemo } from 'react'

const ChartsPanel = () => {
  const mockData = useMemo(
    () => [
      { label: 'Energy (GWh)', value: 1.8 },
      { label: 'Peak Tension (MN)', value: 12.4 },
      { label: 'Glider Throughput (/h)', value: 420 }
    ],
    []
  )

  return (
    <div className="controls-grid">
      {mockData.map((item) => (
        <div key={item.label} className="metric">
          <div>{item.label}</div>
          <strong>{item.value}</strong>
        </div>
      ))}
      <div className="metric">
        Placeholder charts will visualize structural loads, energy usage, and economics over time.
      </div>
    </div>
  )
}

export default ChartsPanel
